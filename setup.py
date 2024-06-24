"""
DPol plugins for SMIRNOFF.
"""
from setuptools import setup

setup(
    name="dpol_plugin",
    packages=["dpol_plugin"],
    version="0.0.0",
    include_package_data=True,
    entry_points={
        "openff.toolkit.plugins.handlers": [
            "DPolMultipoleHandler = dpol_plugin.nonbonded:DPolMultipoleHandler",
            "DPolPolarizabilityHandler = dpol_plugin.nonbonded:DPolPolarizabilityHandler",
        ],
        "openff.interchange.plugins.collections": [
            "DPolCollection = dpol_plugin.nonbonded:DPolCollection",
        ],
    },
)
