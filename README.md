# semp
------------
The Starshade ElectroMagnetic Propagation, or _semp_, is a package that simulates the propagation of light through narrow regions of a starshade. _semp_ relies on the _meep_ software, which uses a Finite Difference Time Domain integration scheme
to solve Maxwell's equations.

Getting started
---------------------
To clone a copy and install: ::

    git clone https://github.com/harnessa/semp.git
    cd semp
    python setup.py install

You'll need to set an environment variable which points to the _semp_ install directory, e.g. (in bash):

    export SEMP=${HOME}/repos/semp

Dependencies
--------------------
You will need:

- `meep <https://meep.readthedocs.io/>`
- `numpy <http://www.numpy.org/>`
- `h5py <http://www.h5py.org>`
- `matplotlib <https://pypi.org/project/matplotlib/>`
- `mpi4py <http://mpi4py.scipy.org/>`
- `pytest <https://pypi.org/project/pytest/>`

Documentation
--------------
No public documentation yet :(

Contributors
------------
Primary author: Anthony Harness (Princeton University)
