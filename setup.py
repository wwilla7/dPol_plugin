"""
DPol plugins for SMIRNOFF.
"""
from setuptools import setup

setup(
    name="mpid_plugin",
    packages=["mpid_plugin"],
    version="0.0.0",
    include_package_data=True,
    entry_points={
        "openff.toolkit.plugins.handlers": [
            "DPolMultipoleHandler = mpid_plugin.nonbonded:DPolMultipoleHandler",
            "DPolPolarizabilityHandler = mpid_plugin.nonbonded:DPolPolarizabilityHandler",
        ],
        "openff.interchange.plugins.collections": [
            "DPolCollection = mpid_plugin.nonbonded:DPolCollection",
        ],
    },
)
