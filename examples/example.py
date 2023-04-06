from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule
from openff.units import unit
from openmm import XmlSerializer

from mpid_plugin.nonbonded import (
    MPIDCollection,
    MPIDMultipoleHandler,
    MPIDPolarizabilityHandler,
)

multipole_handler = MPIDMultipoleHandler(skip_version_check=True)
polarize_handler = MPIDPolarizabilityHandler(skip_version_check=True)

multipol_params = {
    "[#1]-[#8X2H2+0:1]-[#1]": {
        "c0": -0.834 * unit.elementary_charge,
        },
    "[#1:1]-[#8X2H2+0]-[#1]": {
        "c0": 0.417 * unit.elementary_charge,
        }}

polarizability_params = {
    "[#1:1]": {
        "polarizabilityXX": 0.163 * unit.angstrom**3,
        "polarizabilityYY": 0.163 * unit.angstrom**3,
        "polarizabilityZZ": 0.163 * unit.angstrom**3,
        "thole": 8.0* unit.dimensionless}, 
    "[#6:1]": {
        "polarizabilityXX": 1.642 * unit.angstrom**3,
        "polarizabilityYY": 1.642 * unit.angstrom**3,
        "polarizabilityZZ": 1.642 * unit.angstrom**3,
        "thole": 8.0* unit.dimensionless}, 
    "[#7:1]": {
        "polarizabilityXX": 1.042 * unit.angstrom**3,
        "polarizabilityYY": 1.042 * unit.angstrom**3,
        "polarizabilityZZ": 1.042 * unit.angstrom**3,
        "thole": 8.0* unit.dimensionless}, 
    "[#8:1]": {
        "polarizabilityXX": 0.642 * unit.angstrom**3,
        "polarizabilityYY": 0.642 * unit.angstrom**3,
        "polarizabilityZZ": 0.642 * unit.angstrom**3,
        "thole": 8.0* unit.dimensionless}, 
        }


for smirks, values in multipol_params.items():
    multipole_handler.add_parameter(
        {
            "smirks": smirks,
            "c0": values["c0"],
        }
    )

for smirks, values in polarizability_params.items():
    polarize_handler.add_parameter(
        {
            "smirks": smirks,
            "polarizabilityXX": values["polarizabilityXX"],
            "polarizabilityYY": values["polarizabilityYY"],
            "polarizabilityZZ": values["polarizabilityZZ"],
            "thole": values["thole"],
        }
    )

force_field = ForceField("OpenFF-2.0.0.offxml", load_plugins=True)

force_field.register_parameter_handler(multipole_handler)
force_field.register_parameter_handler(polarize_handler)

force_field.to_file("forcefield.offxml")

offmol = Molecule.from_mapped_smiles("[H:2][O:1][H:3]")
offmol.generate_conformers(n_conformers=1)
topology = offmol.to_topology()
position = offmol.conformers[0]

matches = force_field.label_molecules(topology)

print(*matches[0]["MPIDMultipole"].items())
print(*matches[0]["MPIDPolarizability"].items())


ret = Interchange.from_smirnoff(force_field, topology, position)
ret.box = None

system = ret.to_openmm()

with open("system.xml", "w") as file:
    file.write(XmlSerializer.serialize(system))
