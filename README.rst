==========================================
OpenFF SMIRNOFF plugin for polarizability
==========================================

Warning
-------
**This code is experimental and not suitable for production.**

Dependencies
------------
- `openff-toolkit <https://github.com/openforcefield/openff-toolkit>`_
- `openff-interchange >v0.3.0rc1 <https://github.com/openforcefield/openff-interchange>`_
- `mpidplugin <https://github.com/andysim/MPIDOpenMMPlugin>`_
- `openmm >7 <https://github.com/openmm/openmm>`_


Installation
------------

* openff-interchange
    | Because there are two parameter handlers in one collection, we use a slightly modified version of interchange until the new feature gets an official support. `Reference <https://github.com/openforcefield/openff-interchange/pull/648>`_

.. code-block:: sh

    git clone git@github.com:wwilla7/openff-interchange.git
    git checkout pol
    cd openff-interchange
    conda env -f devtools/conda-envs/dev_env.yaml
    pip install -e .

* MPID_Plugin
    | Experimental codes
.. code-block:: sh

    git clone git@github.com:wwilla7/MPID_plugin.git
    cd MPID_plugin
    pip install -e .

* MPIDOpenMMPlugin
    | Follow `instructions <https://github.com/andysim/MPIDOpenMMPlugin>`_
    | Maybe switch to use AMOEBA Force for less dependencies.


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