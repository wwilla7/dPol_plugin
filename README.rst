==========================================
OpenFF SMIRNOFF plugin for polarizability
==========================================

**This code is experimental and not suitable for production.**

Dependencies
------------
- `openff-toolkit <https://github.com/openforcefield/openff-toolkit>`_
- `openff-interchange >v0.3.0rc1 <https://github.com/openforcefield/openff-interchange>`_
- `openmm >=8 <https://github.com/openmm/openmm>`_


Installation
------------

* dPol_Plugin
    | We use the OpenMM Amoeba Force and modified system to compute direct polarization
.. code-block:: sh

    git clone git@github.com:wwilla7/dPol_plugin.git
    cd dPol_plugin
    pip install -e .


Examples
----------

* `example.py <examples/example.py>`_ provides an example to use this plugin.

* `forcefield.xml <examples/forcefield.offxml>`_ is the output openff-style force field file.

* `system.xml <examples/system.xml>`_ is the OpenMM system created by ``openff-interchange``


Copyright
---------

MIT License

Copyright (c) 2023 Open Force Field Initiative

.. include::
    LICENSE