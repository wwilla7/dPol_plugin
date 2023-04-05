"""
MPID plugins for SMIRNOFF.
"""
from setuptools import setup

setup(
    name="mpid_plugin",
    packages=["mpid_plugin"],
    version="0.0.0",
    include_package_data=True,
    entry_points={
        "openff.toolkit.plugins.handlers": [
            "MPIDMultipoleHandler = mpid_plugin.nonbonded:MPIDMultipoleHandler",
            "MPIDPolarizabilityHandler = mpid_plugin.nonbonded:MPIDPolarizabilityHandler",
        ],
        "openff.interchange.plugins.collections": [
            "MPIDCollection = mpid_plugin.nonbonded:MPIDCollection",
        ],
    },
)
