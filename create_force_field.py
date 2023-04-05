from openff.toolkit import ForceField, Molecule
from openff.units import unit

from mpid_plugin.nonbonded import (
    MPIDCollection,
    MPIDMultipoleHandler,
    MPIDPolarizabilityHandler,
)

multipole_handler = MPIDMultipoleHandler(skip_version_check=True)
polarize_handler = MPIDPolarizabilityHandler(skip_version_check=True)

parameters = {
    "[#1]-[#8X2H2+0:1]-[#1]": {
        "c0": -0.834 * unit.elementary_charge,
        "polarizabilityXX": 0.00088 * unit.nanometer**3,
        "polarizabilityYY": 0.00088 * unit.nanometer**3,
        "polarizabilityZZ": 0.00088 * unit.nanometer**3,
        "thole": 8.0 * unit.dimensionless,
    },
    "[#1:1]-[#8X2H2+0]-[#1]": {
        "c0": 0.417 * unit.elementary_charge,
        # hydrogen has no polarizability, these should just not exist
        "polarizabilityXX": 0.00088 * unit.nanometer**3,
        "polarizabilityYY": 0.00088 * unit.nanometer**3,
        "polarizabilityZZ": 0.00088 * unit.nanometer**3,
        "thole": 0.0 * unit.dimensionless,
    },
}

for smirks, values in parameters.items():
    multipole_handler.add_parameter(
        {
            "smirks": smirks,
            "c0": values["c0"],
        }
    )
    polarize_handler.add_parameter(
        {
            "smirks": smirks,
            "polarizabilityXX": values["polarizabilityXX"],
            "polarizabilityYY": values["polarizabilityYY"],
            "polarizabilityZZ": values["polarizabilityZZ"],
            "thole": values["thole"],
        }
    )

force_field = ForceField(load_plugins=True)

force_field.register_parameter_handler(multipole_handler)
force_field.register_parameter_handler(polarize_handler)

force_field.to_file("mpid.offxml")


topology = Molecule.from_mapped_smiles("[H:2][O:1][H:3]").to_topology()

matches = force_field.label_molecules(topology)

print(*matches[0]["MPIDMultipole"].items())
print(*matches[0]["MPIDPolarizability"].items())
