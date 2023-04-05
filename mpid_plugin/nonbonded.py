from typing import TYPE_CHECKING, Dict, List, Literal, Optional, Set, Tuple, Union

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

if TYPE_CHECKING:
    from openff.interchange.models import VirtualSiteKey

# https://docs.openforcefield.org/projects/interchange/en/v0.3.0-staging/using/plugins.html


class MultipoleKey(TopologyKey):
    pass


class PolarizabilityKey(TopologyKey):
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
    type: Literal["MPID"] = "MPID"

    expression: str = ""

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

    def store_potentials(
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
        polarizability_handler = polarizability_handler.find_matches(topology)

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
                associated_handler=parameter_handler.TAGNAME,
            )
            self.key_map[topology_key] = potential_key

        for key, val in polarizability_handler.items():
            topology_key = PolarizabilityKey(atom_indices=key)
            potential_key = PotentialKey(
                id=val.parameter_type.smirks,
                associated_handler=parameter_handler.TAGNAME,
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

        for potential_key in self.key_map.values():
            if potential_key.associated_handler == "MPIDMultipole":
                smirks = potential_key.id
                parameter = multipole_handler.parameters[smirks]

                self.potentials[potential_key] = Potential(
                    parameters={"c0": parameter.c0}
                )

            if potential_key.associated_handler == "MPIDPolarizability":
                smirks = potential_key.id
                parameter = polarizability_handler.parameters[smirks]

                self.potentials[potential_key] = Potential(
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
        parameter_handler: List[ParameterAttribute],
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
            nonbonded_force.setParticleParameters(
                particle_index, 0.0 * unit.elementary_charge, sigma, epsilon
            )

        # Pesudocode from here to end of file!
        # Create the MPID force
        mpid_collection = interchange.collections["MPID"]

        mpid_force = openmm.MPIDForce(coulomb14scale=mpid_collection.coulomb14scale)
        system.addForce(mpid_force)

        # Set the multipole and polarizability parameters on the force
        for topology_key, potential_key in mpid_collection.potentials.items():
            openff_particle_index = topology_key.atom_indices[0]
            openmm_particle_index = particle_map[openff_particle_index]

            if isinstance(topology_key, MultipoleKey):
                multipole: Potential = mpid_collection.potentials[potential_key]
                # Set multipole on multipole force using OpenMM particle index,
                # c0 from `potential_key.parameters['c0']`
                mpid_force.addMultipole(...)

            if isinstance(topology_key, PolarizabilityKey):
                polarizability: Potential = mpid_collection.potentials[potential_key]
                # Set polarizability on multipole force using OpenMM particle index,
                # polarizabilityXX, polarizabilityYY, polarizabilityZZ, thole from
                # `potential_key.parameters['polarizabilityXX']`, etc.
                mpid_force.addPolarizability(...)

        # Plus whatever other housekeeping needs to happen with the OpenMM forces
