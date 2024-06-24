from mpid_plugin.nonbonded import (DPolCollection, DPolMultipoleHandler,
                                   DPolPolarizabilityHandler)
from openeye import oechem
from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule, Topology
from openff.toolkit.utils.exceptions import DuplicateParameterError
from openff.units import unit
from openmm import XmlSerializer

multipole_handler = DPolMultipoleHandler(skip_version_check=True)
polarize_handler = DPolPolarizabilityHandler(skip_version_check=True)

water_parameters = {
    "Bonds": {
        "[#1:1]-[#8X2H2+0:2]-[#1]": {
            "length": 0.1037789311313845 * unit.nanometer,
            "k": 462750.4 * unit.kilojoule_per_mole / unit.nanometer / unit.nanometer,
        }
    },
    "Angles": {
        "[*:1]-[#8:2]-[*:3]": {
            "angle": 1.9527272143303327 * unit.radian,
            "k": 836.8 * unit.kilojoule_per_mole / unit.radian / unit.radian,
        }
    },
    "vdW": {
        "[#1]-[#8X2H2+0:1]-[#1]": {
            "sigma": 0.31897404781098154 * unit.nanometer,
            "epsilon": 0.6389151692355459 * unit.kilojoule_per_mole,
        },
        "[#1:1]-[#8X2H2+0]-[#1]": {
            "sigma": 1.0 * unit.nanometer,
            "epsilon": 0.0 * unit.kilojoule_per_mole,
        },
    },
    "Electrostatics": {
        "[#1]-[#8X2H2+0:1]-[#1]": {"c0": -0.6870676714866755 * unit.elementary_charge},
        "[#1:1]-[#8X2H2+0]-[#1]": {"c0": 0.34353383574333773 * unit.elementary_charge},
    },
    "Polarizability": {
        "[#1:1]-[#8X2H2+0]-[#1]": {
            "polarizability": 0.00012531035951495126 * unit.nanometer**3
        },
        "[#1]-[#8X2H2+0:1]-[#1]": {
            "polarizability": 0.0008987586822022096 * unit.nanometer**3
        },
    },
}

multipol_params = water_parameters["Electrostatics"]

polarizability_params = {
    "[#1:1]": {
        "polarizability": 0.163 * unit.angstrom**3,
    },
    "[#6:1]": {
        "polarizability": 1.642 * unit.angstrom**3,
    },
    "[#7:1]": {
        "polarizability": 1.042 * unit.angstrom**3,
    },
    "[#8:1]": {
        "polarizability": 0.642 * unit.angstrom**3,
    },
} | water_parameters["Polarizability"]


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
            "polarizability": values["polarizability"],
        }
    )

force_field = ForceField("OpenFF-2.2.0.offxml", load_plugins=True)


force_field.register_parameter_handler(multipole_handler)
force_field.register_parameter_handler(polarize_handler)
force_field.deregister_parameter_handler("Constraints")

for p in ["Bonds", "Angles", "vdW"]:
    for s, v in water_parameters[p].items():
        try:
            force_field[p].add_parameter({"smirks": s} | v)
        except DuplicateParameterError:
            old = force_field[p].get_parameter({"smirks": s})[0]
            force_field[p].parameters.remove(old)
            force_field[p].add_parameter({"smirks": s} | v)

force_field.to_file("forcefield.offxml")

molecule = Molecule.from_smiles("O")
molecule.generate_conformers(n_conformers=1)
topology = molecule.to_topology()
position = molecule.conformers[0]

ret = Interchange.from_smirnoff(
    force_field=force_field, topology=topology, positions=position
)
ret.to_openmm()

system = ret.to_openmm()

with open("system.xml", "w") as file:
    file.write(XmlSerializer.serialize(system))
