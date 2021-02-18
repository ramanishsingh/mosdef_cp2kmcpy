Installation
============
Conda installation is recommended. The package can be downloaded from `here <https://github.com/ramanishsingh/cp2kmdpy>`_.


Conda environment
-----------------
cp2kmdpy has been built and tested on python 3.7.
You can use the following commands to create a virtual conda environment with all the cp2kmdpy dependencies.


>>> conda create -n mosdef36 -c mosdef -c conda-forge -c omnia python=3.6 mbuild foyer signac signac-flow pandas
>>> conda activate mosdef36
>>> conda install -c conda-forge openbabel
>>> pip install unyt
>>> pip install ele
>>> git clone https://github.com/ramanishsingh/cp2kmdpy.git
>>> cd cp2kmdpy
>>> pip install -e .
>>> git clone https://github.com/ramanishsingh/mosdef_cp2k_writer
>>> cd mosdef_cp2k_writer
>>> pip install -e .



