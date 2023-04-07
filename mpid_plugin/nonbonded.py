from typing import TYPE_CHECKING, Dict, List, Literal, Optional, Set, Tuple, Union
from pydantic import Field

import openmm
from openff.interchange import Interchange
from openff.interchange.components.potentials import Potential
from openff.interchange.models import PotentialKey, TopologyKey
from openff.interchange.smirnoff._nonbonded import SMIRNOFFCollection
from openff.models.types import FloatQuantity
from openff.toolkit import Topology
from openff.toolkit.typing.engines.smirnoff.parameters import (
    ParameterAttribute,
    ParameterHandler,
    ParameterType,
    _allow_only,
    unit,
)

import mpidplugin
from mpidplugin import MPIDForce
from openmm.app import ForceField
from openmm.app import NoCutoff, PME, LJPME
NONBONDED_METHODS = { 0 : NoCutoff, 4: PME, 5:LJPME}
from openmm.app.forcefield import AmoebaVdwGenerator
import numpy as np

if TYPE_CHECKING:
    from openff.interchange.models import VirtualSiteKey

# https://docs.openforcefield.org/projects/interchange/en/v0.3.0-staging/using/plugins.html


class MultipoleKey(TopologyKey):
    name = Literal["Multipole"]
    pass


class PolarizabilityKey(TopologyKey):
    name = Literal["Polarizability"]
    pass


class MPIDMultipoleHandler(ParameterHandler):
    """
    <MPIDForce coulomb14scale="1.0" >
      <Multipole type="OW" c0="-0.834" />
      <Multipole type="HW" c0="0.417" />
      <Multipole type="HA3" c0="0.09" />
      <Multipole type="CT3" c0="-0.27" />
      <Polarize type="OW" polarizabilityXX="0.00088" polarizabilityYY="0.00088" polarizabilityZZ="0.00088" thole="8.0"/>
      <Polarize type="CT3" polarizabilityXX="0.00068" polarizabilityYY="0.00068" polarizabilityZZ="0.00068" thole="8.0"/>
    </MPIDForce>
    """

    class MPIDMultipole(ParameterType):
        _VALENCE_TYPE = "Atom"
        _ELEMENT_NAME = "Multipole"

        # This is just the charge, which is handled here instead of the NonbondedForce
        c0 = ParameterAttribute(default=None, unit=unit.elementary_charge)

    _TAGNAME = "MPIDMultipole"
    _INFOTYPE = MPIDMultipole

    # https://github.com/andysim/MPIDOpenMMPlugin/blob/43450d73e567772e8892cabf9dde7f6c34913e4e/examples/ethane_water_charge_only/ethane_water.xml#L57
    coulomb14scale = ParameterAttribute(default=1.0, converter=float)


class MPIDPolarizabilityHandler(ParameterHandler):
    class MPIDPolarizability(ParameterType):
        _VALENCE_TYPE = "Atom"
        _ELEMENT_NAME = "Polarizability"

        polarizabilityXX = ParameterAttribute(default=None, unit=unit.nanometer**3)
        polarizabilityYY = ParameterAttribute(default=None, unit=unit.nanometer**3)
        polarizabilityZZ = ParameterAttribute(default=None, unit=unit.nanometer**3)

        thole = ParameterAttribute(default=None, unit=unit.dimensionless)

    _TAGNAME = "MPIDPolarizability"
    _INFOTYPE = MPIDPolarizability

    # https://github.com/andysim/MPIDOpenMMPlugin/blob/43450d73e567772e8892cabf9dde7f6c34913e4e/examples/ethane_water_charge_only/ethane_water.xml#L57
    coulomb14scale = ParameterAttribute(default=1.0, converter=float)


class MPIDCollection(SMIRNOFFCollection):
    is_plugin = True

    type: Literal["MPID"] = "MPID"

    expression: str = "Direct Polarization"

    coulomb14scale: float = Field(
        1.0, description="The scaling factor applied to 1-4 interactions"
    )

    @classmethod
    def allowed_parameter_handlers(cls):
        return [MPIDMultipoleHandler, MPIDPolarizabilityHandler]

    @classmethod
    def supported_parameters(cls):
        return (
            "smirks",
            "id",
            "c0",
            "polarizabilityXX",
            "polarizabilityYY",
            "polarizabilityZZ",
            "thole",
        )

    def store_matches(
        self,
        parameter_handler: List[ParameterHandler],
        topology: "Topology",
    ) -> None:
        """Populate self.key_map with key-val pairs of [TopologyKey, PotentialKey]."""
        self.key_map: Dict[TopologyKey, PotentialKey] = dict()

        # Assume there are two parameter handlers, one for multipole and one for polarizability
        multipole_handler = [
            x for x in parameter_handler if x.TAGNAME == "MPIDMultipole"
        ][0]
        polarizability_handler = [
            x for x in parameter_handler if x.TAGNAME == "MPIDPolarizability"
        ][0]

        # The multipole stores charges, so assume all atoms have a multipole
        multipole_matches = multipole_handler.find_matches(topology)

        # not all atoms have polarizbility
        polarizability_matches = polarizability_handler.find_matches(topology)

        # WW has custom charges (stored as multipole parameters) using a custom model,
        # which has potentially many smirks-to-charge mappings. By contrast, there are
        # only polarizability terms on CHON atoms, and (for now) all C have the same
        # polarizability parameters, all hydrogens, etc. (In the future there may be more,
        # or more complex SMIRKS patterns.) So need to store the different types of keys,
        # cannot collapse them.
        for key, val in multipole_matches.items():
            topology_key = MultipoleKey(atom_indices=key)
            potential_key = PotentialKey(
                id=val.parameter_type.smirks,
                associated_handler=multipole_handler.TAGNAME,
            )
            self.key_map[topology_key] = potential_key

        for key, val in polarizability_matches.items():
            topology_key = PolarizabilityKey(atom_indices=key)
            potential_key = PotentialKey(
                id=val.parameter_type.smirks,
                associated_handler=polarizability_handler.TAGNAME,
            )
            self.key_map[topology_key] = potential_key

    def store_potentials(self, parameter_handler: List[ParameterHandler]):
        multipole_handler = [
            x for x in parameter_handler if x.TAGNAME == "MPIDMultipole"
        ][0]
        polarizability_handler = [
            x for x in parameter_handler if x.TAGNAME == "MPIDPolarizability"
        ][0]

        self.coulomb14scale = parameter_handler[0].coulomb14scale

        for topology_key, potential_key in self.key_map.items():
            if potential_key.associated_handler == "MPIDMultipole":
                smirks = potential_key.id
                parameter = multipole_handler.parameters[smirks]

                self.potentials[topology_key] = Potential(
                    parameters={"c0": parameter.c0}
                )

            if potential_key.associated_handler == "MPIDPolarizability":
                smirks = potential_key.id
                parameter = polarizability_handler.parameters[smirks]

                self.potentials[topology_key] = Potential(
                    parameters={
                        "polarizabilityXX": parameter.polarizabilityXX,
                        "polarizabilityYY": parameter.polarizabilityYY,
                        "polarizabilityZZ": parameter.polarizabilityZZ,
                        "thole": parameter.thole,
                    }
                )

    @classmethod
    def create(
        cls,
        parameter_handler: List[ParameterHandler],
        topology: Topology,
    ):
        # Assume the two handlers have the same coluomb14scale
        handler = cls(
            coulomb14scale=parameter_handler[0].coulomb14scale,
        )

        handler.store_matches(parameter_handler=parameter_handler, topology=topology)
        handler.store_potentials(parameter_handler=parameter_handler)

        return handler

    @classmethod
    def check_openmm_requirements(cls, combine_nonbonded_forces: bool) -> None:
        """Later, when setting charges to 0, we will assume that there's just one NonbondedForce."""
        assert combine_nonbonded_forces

    def modify_openmm_forces(
        self,
        interchange: Interchange,
        system: openmm.System,
        add_constrained_forces: bool,
        constrained_pairs: Set[Tuple[int, ...]],
        particle_map: Dict[Union[int, "VirtualSiteKey"], int],
    ):
        # Set the charges on the nonbonded force to be zero
        nonbonded_force = [
            force
            for force in system.getForces()
            if isinstance(force, openmm.NonbondedForce)
        ][0]

        for particle_index in range(nonbonded_force.getNumParticles()):
            _, sigma, epsilon = nonbonded_force.getParticleParameters(particle_index)
            nonbonded_force.setParticleParameters(particle_index, 0.0, sigma, epsilon)
        
        nonbonded_method = NONBONDED_METHODS[nonbonded_force.getNonbondedMethod()]
        cutoff_distance = nonbonded_force.getCutoffDistance()

        # Pesudocode from here to end of file!
        # First attempt to create MPIDForce
        # Create the MPID force
        methodMap = {NoCutoff:MPIDForce.NoCutoff,
                     PME:MPIDForce.PME,
                     LJPME:MPIDForce.PME}
        mpid_collection = interchange.collections["MPID"]

        mpid_force = MPIDForce()
        
        mpid_force.setNonbondedMethod(methodMap[nonbonded_method])

        mpid_force.setCutoffDistance(cutoff_distance)
        mpid_force.setPolarizationType(MPIDForce.Direct)
        mpid_force.set14ScaleFactor(mpid_collection.coulomb14scale)

        # Every atom has partial charge but not every atom has polarizability parameter
        # Let's start with all of them have both charges and polarizabilities
        # No multipole parameters

        n_particles = system.getNumParticles()

        # Because we use direct polarization and monopole, we can preset these parameters
        parameter_maps = {
            i: {
                "dipole": np.zeros(3).tolist(),
                "quadrupole": np.zeros(6).tolist(),
                "octopole": np.zeros(10).tolist(),
                "axisType": 5,
                "thole": 8.0,
            }
            for i in range(n_particles)
        }

        # Set the multipole and polarizability parameters on the force
        for topology_key, potential_key in mpid_collection.potentials.items():
            openff_particle_index = topology_key.atom_indices[0]
            openmm_particle_index = particle_map[openff_particle_index]

            if isinstance(topology_key, MultipoleKey):
                # multipole: Potential = mpid_collection.potentials[potential_key]
                # Set multipole on multipole force using OpenMM particle index,
                # c0 from `potential_key.parameters['c0']`
                parameter_maps[openmm_particle_index][
                    "charge"
                ] = potential_key.parameters["c0"].m_as(unit.elementary_charge)

            if isinstance(topology_key, PolarizabilityKey):
                # polarizability: Potential = mpid_collection.potentials[potential_key]
                # Set polarizability on multipole force using OpenMM particle index,
                # polarizabilityXX, polarizabilityYY, polarizabilityZZ, thole from
                # `potential_key.parameters['polarizabilityXX']`, etc.
                parameter_maps[openmm_particle_index]["polarizability"] = [
                    potential_key.parameters["polarizabilityXX"].m_as(
                        unit.nanometer**3
                    ),
                    potential_key.parameters["polarizabilityYY"].m_as(
                        unit.nanometer**3
                    ),
                    potential_key.parameters["polarizabilityZZ"].m_as(
                        unit.nanometer**3
                    ),
                ]

                parameter_maps[openmm_particle_index]["thole"] = float(
                    potential_key.parameters["thole"].magnitude
                )

        ## Add multipoles
        ## Find covalentMap
        omm_ff = ForceField()
        data = omm_ff._SystemData(interchange.topology.to_openmm())
        bonded12ParticleSets = AmoebaVdwGenerator.getBondedParticleSets(system, data)

        bonded13ParticleSets = []
        for i in range(len(data.atoms)):
            bonded13Set = set()
            bonded12ParticleSet = bonded12ParticleSets[i]
            for j in bonded12ParticleSet:
                bonded13Set = bonded13Set.union(bonded12ParticleSets[j])

        # remove 1-2 and self from set

            bonded13Set = bonded13Set - bonded12ParticleSet
            selfSet = set()
            selfSet.add(i)
            bonded13Set = bonded13Set - selfSet
            bonded13Set = set(sorted(bonded13Set))
            bonded13ParticleSets.append(bonded13Set)

        # 1-4

        bonded14ParticleSets = []
        for i in range(len(data.atoms)):
            bonded14Set = set()
            bonded13ParticleSet = bonded13ParticleSets[i]
            for j in bonded13ParticleSet:
                bonded14Set = bonded14Set.union(bonded12ParticleSets[j])

        # remove 1-3, 1-2 and self from set

            bonded14Set = bonded14Set - bonded12ParticleSets[i]
            bonded14Set = bonded14Set - bonded13ParticleSet
            selfSet = set()
            selfSet.add(i)
            bonded14Set = bonded14Set - selfSet
            bonded14Set = set(sorted(bonded14Set))
            bonded14ParticleSets.append(bonded14Set)

        for particle in range(n_particles):
            mpid_force.addMultipole(
                parameter_maps[particle]["charge"],
                parameter_maps[particle]["dipole"],
                parameter_maps[particle]["quadrupole"],
                parameter_maps[particle]["quadrupole"],
                parameter_maps[particle]["axisType"],
                -1,
                -1,
                -1,
                parameter_maps[particle]["thole"],
                parameter_maps[particle]["polarizability"],
            )
        ## TODO: Add `CovalentMap`
            mpid_force.setCovalentMap(particle, MPIDForce.Covalent12, tuple(bonded12ParticleSets[particle]))
            mpid_force.setCovalentMap(particle, MPIDForce.Covalent13, tuple(bonded13ParticleSets[particle]))
            mpid_force.setCovalentMap(particle, MPIDForce.Covalent14, tuple(bonded14ParticleSets[particle]))

        system.addForce(mpid_force)


        # Plus whatever other housekeeping needs to happen with the OpenMM forces
