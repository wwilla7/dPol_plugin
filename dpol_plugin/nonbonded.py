from typing import TYPE_CHECKING, Dict, List, Literal, Set, Tuple, Union

import numpy as np
import openmm
from openff.interchange import Interchange
from openff.interchange.components.potentials import Potential
from openff.interchange.models import (
    PotentialKey,
    TopologyKey,
    LibraryChargeTopologyKey,
)
from openff.interchange.smirnoff._nonbonded import SMIRNOFFCollection
from openff.toolkit import Topology
from openff.toolkit.typing.engines.smirnoff.parameters import (
    ParameterAttribute,
    ParameterHandler,
    ParameterType,
    unit,
)
from openmm import AmoebaMultipoleForce
from openmm.app import LJPME, PME, ForceField, NoCutoff
from openmm.app.forcefield import AmoebaVdwGenerator

NONBONDED_METHODS = {0: NoCutoff, 4: PME}

if TYPE_CHECKING:
    from openff.interchange.models import VirtualSiteKey

# https://docs.openforcefield.org/projects/interchange/en/v0.3.0-staging/using/plugins.html


class MultipoleKey(LibraryChargeTopologyKey):
    # name = Literal["Multipole"]
    type: Literal["Multipole"] = "Multipole"


class PolarizabilityKey(TopologyKey):
    # name = Literal["Polarizability"]
    type: Literal["Polarizability"] = "Polarizability"


class DPolMultipoleHandler(ParameterHandler):
    """
    <AmoebaMultipoleForce direct11Scale="0.0" direct12Scale="1.0" direct13Scale="1.0" direct14Scale="1.0" mpole12Scale="0.0"
    mpole13Scale="0.0" mpole14Scale="0.5" mpole15Scale="1.0" mutual11Scale="1.0" mutual12Scale="1.0" mutual13Scale="1.0"
    mutual14Scale="1.0" polar12Scale="0.0" polar13Scale="0.0" polar14Intra="0.5" polar14Scale="1.0" polar15Scale="1.0">
      <Multipole type="301" kz="0" c0="-0.6870676714866755"
                                   d1="0.0" d2="0.0" d3="0.0" q11="0.0" q21="0.0" q22="0.0" q31="0.0" q32="0.0" q33="0.0"/>
      <Multipole type="302" kz="0" c0="0.34353383574333773"
                                   d1="0.0" d2="0.0" d3="0.0" q11="0.0" q21="0.0" q22="0.0" q31="0.0" q32="0.0" q33="0.0"/>
      <Polarize type="301" polarizability="0.0008987586822022096" thole="0.0" pgrp1="302"/>
      <Polarize type="302" polarizability="0.00012531035951495126" thole="0.0" pgrp1="301"/>
    </AmoebaMultipoleForce>
    """

    class DPolMultipole(ParameterType):
        _VALENCE_TYPE = "Atom"
        _ELEMENT_NAME = "Multipole"

        # This is just the charge, which is handled here instead of the NonbondedForce
        c0 = ParameterAttribute(default=None, unit=unit.elementary_charge)

    _TAGNAME = "DPolMultipole"
    _INFOTYPE = DPolMultipole


class DPolPolarizabilityHandler(ParameterHandler):
    class DPolPolarizability(ParameterType):
        _VALENCE_TYPE = "Atom"
        _ELEMENT_NAME = "Polarizability"

        polarizability = ParameterAttribute(default=None, unit=unit.nanometer**3)

        # direct polarization model doesn't need thole damping, set to zero
        # thole = ParameterAttribute(default=0.0, unit=unit.dimensionless)

    _TAGNAME = "DPolPolarizability"
    _INFOTYPE = DPolPolarizability


class DPolCollection(SMIRNOFFCollection):
    is_plugin: bool = True
    acts_as: str = "DPol"

    type: Literal["DPol"] = "DPol"

    expression: str = "Direct Polarization"

    @classmethod
    def allowed_parameter_handlers(cls):
        return [DPolMultipoleHandler, DPolPolarizabilityHandler]

    @classmethod
    def supported_parameters(cls):
        return (
            "smirks",
            "id",
            "c0",
            "polarizability",
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
            x for x in parameter_handler if x.TAGNAME == "DPolMultipole"
        ][0]
        polarizability_handler = [
            x for x in parameter_handler if x.TAGNAME == "DPolPolarizability"
        ][0]

        # The multipole stores charges, so assume all atoms have a multipole
        multipole_matches = multipole_handler.find_matches(topology)

        polarizability_matches = polarizability_handler.find_matches(topology)

        # WW has custom charges (stored as multipole parameters) using a custom model,
        # which has potentially many smirks-to-charge mappings. By contrast, there are
        # only polarizability terms on CHON atoms, and (for now) all C have the same
        # polarizability parameters, all hydrogens, etc. (In the future there may be more,
        # or more complex SMIRKS patterns.) So need to store the different types of keys,
        # cannot collapse them.
        for key, val in multipole_matches.items():
            topology_key = MultipoleKey(this_atom_index=key[0])
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
            x for x in parameter_handler if x.TAGNAME == "DPolMultipole"
        ][0]
        polarizability_handler = [
            x for x in parameter_handler if x.TAGNAME == "DPolPolarizability"
        ][0]

        for topology_key, potential_key in self.key_map.items():
            if potential_key.associated_handler == "DPolMultipole":
                smirks = potential_key.id
                parameter = multipole_handler.parameters[smirks]

                self.potentials[topology_key] = Potential(
                    parameters={"c0": parameter.c0}
                )

            if potential_key.associated_handler == "DPolPolarizability":
                smirks = potential_key.id
                parameter = polarizability_handler.parameters[smirks]

                self.potentials[topology_key] = Potential(
                    parameters={
                        "polarizability": parameter.polarizability,
                    }
                )

    @classmethod
    def create(
        cls,
        parameter_handler: List[ParameterHandler],
        topology: Topology,
    ):
        handler = cls()
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
        # First attempt to create DPolForce
        # Create the DPol force
        methodMap = {
            NoCutoff: AmoebaMultipoleForce.NoCutoff,
            PME: AmoebaMultipoleForce.PME,
        }
        mpid_collection = interchange.collections["DPol"]

        dpol_force = AmoebaMultipoleForce()

        dpol_force.setNonbondedMethod(methodMap[nonbonded_method])

        dpol_force.setCutoffDistance(cutoff_distance)
        dpol_force.setPolarizationType(AmoebaMultipoleForce.Direct)

        # Every atom has partial charge but not every atom has polarizability parameter
        # Let's start with all of them have both charges and polarizabilities
        # No multipole parameters

        n_particles = system.getNumParticles()

        # Because we use direct polarization and monopole, we can preset these parameters
        parameter_maps = {
            i: {
                "dipole": np.zeros(3).tolist(),
                "quadrupole": np.zeros(9).tolist(),
                # notes on axisType
                # https://github.com/openmm/openmm/blob/71b2a93e09d9d5f110310448c9ff503eb2e71d55/plugins/amoeba/openmmapi/include/openmm/AmoebaMultipoleForce.h#L96
                "axisType": 5,
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
                parameter_maps[openmm_particle_index][
                    "polarizability"
                ] = potential_key.parameters["polarizability"].m_as(unit.nanometer**3)

        # Find covalentMap
        omm_ff = ForceField()
        data = omm_ff._SystemData(interchange.topology.to_openmm())

        ## Reference
        ## https://github.com/openmm/openmm/blob/master/plugins/amoeba/openmmapi/src/AmoebaMultipoleForce.cpp
        ## https://github.com/openmm/openmm/blob/71b2a93e09d9d5f110310448c9ff503eb2e71d55/wrappers/python/openmm/app/forcefield.py#L4909
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

        # 1-5

        bonded15ParticleSets = []
        for i in range(len(data.atoms)):
            bonded15Set = set()
            bonded14ParticleSet = bonded14ParticleSets[i]
            for j in bonded14ParticleSet:
                bonded15Set = bonded15Set.union(bonded12ParticleSets[j])

            # remove 1-4, 1-3, 1-2 and self from set

            bonded15Set = bonded15Set - bonded12ParticleSets[i]
            bonded15Set = bonded15Set - bonded13ParticleSets[i]
            bonded15Set = bonded15Set - bonded14ParticleSet
            selfSet = set()
            selfSet.add(i)
            bonded15Set = bonded15Set - selfSet
            bonded15Set = set(sorted(bonded15Set))
            bonded15ParticleSets.append(bonded15Set)

        for particle in range(n_particles):
            dpol_force.addMultipole(
                charge=parameter_maps[particle]["charge"],
                molecularDipole=parameter_maps[particle]["dipole"],
                molecularQuadrupole=parameter_maps[particle]["quadrupole"],
                axisType=parameter_maps[particle]["axisType"],
                # zaxis, xaxis, yaxis
                multipoleAtomZ=-1,
                multipoleAtomX=-1,
                multipoleAtomY=-1,
                thole=0,
                dampingFactor=0,
                polarity=parameter_maps[particle]["polarizability"],
            )

            dpol_force.setCovalentMap(
                particle,
                AmoebaMultipoleForce.Covalent12,
                tuple(bonded12ParticleSets[particle]),
            )
            dpol_force.setCovalentMap(
                particle,
                AmoebaMultipoleForce.Covalent13,
                tuple(bonded13ParticleSets[particle]),
            )
            dpol_force.setCovalentMap(
                particle,
                AmoebaMultipoleForce.Covalent14,
                tuple(bonded14ParticleSets[particle]),
            )
            dpol_force.setCovalentMap(
                particle,
                AmoebaMultipoleForce.Covalent15,
                tuple(bonded15ParticleSets[particle]),
            )

            # add artificial covalent 15 on nonbonded force
            covalent15 = tuple(bonded15ParticleSets[particle])
            if len(covalent15) > 0:
                for p in covalent15:
                    if p < particle:
                        q1, s1, e1 = nonbonded_force.getParticleParameters(particle)
                        q2, s2, e2 = nonbonded_force.getParticleParameters(p)
                        nonbonded_force.addException(
                            particle1=particle,
                            particle2=p,
                            chargeProd=0 * openmm.unit.elementary_charge**2,
                            sigma=0.5 * (s1 + s2),
                            epsilon=np.sqrt(e1 * e2),
                        )

            # put 1-2, 1-3 in the same polarization group
            polar11 = (
                tuple(bonded12ParticleSets[particle])
                + tuple(bonded13ParticleSets[particle])
                + (particle,)
            )
            dpol_force.setCovalentMap(
                particle, AmoebaMultipoleForce.PolarizationCovalent11, polar11
            )

            polar12 = tuple(bonded14ParticleSets[particle])
            dpol_force.setCovalentMap(
                particle, AmoebaMultipoleForce.PolarizationCovalent12, polar12
            )

            polar13 = tuple(bonded15ParticleSets[particle])
            dpol_force.setCovalentMap(
                particle, AmoebaMultipoleForce.PolarizationCovalent13, polar13
            )

        system.addForce(dpol_force)
