Running Monte Carlo
===========================

This functionality of the mosdef_cp2kmcpy package can be used to run Monte Carlo simulations in CP2K. 

.. autoclass:: cp2kmdpy.md.MD
    :members:

Instantiating ``MD`` first requires the specification of the molecule(s), functional, box , the basis set, and the number of molecules of each type.
It can be done as follows:

.. code-block:: python

  import cp2kmdpy
  md_instance = cp2kmdpy.md.MD(
      molecules=[molecule],
      functional=functional,
      box=box,
      basis_set=basis_set,
      n_molecules=n_molecules
  )

MD class
****************************

molecules
~~~~~~~~~~
``molecules`` is the list of molecules present in the simulation box. Each element of the list must belong to the ``mbuild.Compound`` class.

For example, for a mixture of ethane and propane, ``molecules`` can be defined as follows:

.. code-block:: python

  import mbuild

  ethane = mbuild.load("CC", smiles=True)
  propane = mbuild.load("CCC", smiles=True)
  molecules= [ethane, propane]

functional
~~~~~~~~~~~~
``functional`` is the name of the DFT XC functional to be used for the MD simulation.
``functional`` should be a string.

.. code-block:: python

  functional = 'PBE'


box
~~~~
``box`` is the simulation box which will be populated with the molecules.
``box`` should be an instance of the  ``mbuild.Box`` class.

.. code-block:: python

  import mbuild
  box = mbuild.Box(lengths=[1.0, 1.0, 1.0], angles=[90., 90., 90.])

cutoff
~~~~~~
While using the GPW method in CP2K, the plane wave cutoff needs to be defined.
``cutoff`` is the plane wave cutoff in Ry for optimization simulation. It should be a float.


.. code-block:: python

  cutoff=600

scf_tolerance
~~~~~~~~~~~~~~~
``scf_tolerance`` controls the tolerance of each SCF cycle in the optimization. Should be a float.

.. code-block:: python

  scf_tolerance=1e-7

basis_set
~~~~~~~~~
With ``basis_set``, the basis set to be used for each element present in the simulation can be defined. It is a dictionary with keys being the element symbol and values being the basis set for that particular element.
For example, if the molecule contains carbon and oxygen atoms and DZVP-MOLOPT-SR-GTH basis set is desired, then the ``basis_set`` can be defined as:

.. code-block:: python

  basis_set={'C':'DZVP-MOLOPT-SR-GTH','O':'DZVP-MOLOPT-GTH'}

basis_set_filename
~~~~~~~~~~~~~~~~~~~
The name of the file to look into for the basis set specified. Should be a string. Default value is ``BASIS_MOLOPT``.

potential_filename
~~~~~~~~~~~~~~~~~~~
The name of the file to look into for pseudopotential information. Should be a string. Default value is ``GTH_POTENTIALS``.

fixed_list
~~~~~~~~~~~
During molecular dynamics simulations, sometimes constraining the movement of some atoms is desired. This can be achieved by specifying ``fixed_atoms``.
For example, if atom number 1 to 80 needs to be fixed, it can be specified as follows

.. code-block:: python

  fixed_list='1..80'

periodicity
~~~~~~~~~~~
This attribute controls the periodicity of the box. Default value is ``'XYZ'``
Consider a system which is periodic in only X and Y direction. It can be specified as follows:

.. code-block:: python

  periodicity='XY'

simulation_time
~~~~~~~~~~~~~~~~
Time for which simulation should be run. It should be a ``unyt_quantity``.
For example if the simulation should be run for 1 ps , ``simulation_time`` can be set as:

.. code-block:: python

  import unyt as u

  simulation_time=1*u.ps

The default value of ``simulation_time`` is 10 ps.

time_step
~~~~~~~~~~
Time step for the molecula dynamics simulation. Should be a ``unyt_quantity``.

.. code-block:: python

  import unyt as u

  time_step=1*u.fs

ensemble
~~~~~~~~~~
Ensemble for the molecular dynamics simulation. Must be a string of capital letters.
Deafult value is ``'NVT'``

project_name
~~~~~~~~~~~~~~~
Name of the project. Should be string. Default value is ``'sample_project'``

temperature
~~~~~~~~~~~~
Temperature of the ensemble, if required. Must be specified as a ``unyt_quantity``.

.. code-block:: python

  import unyt as u

  temperature=298*u.K

pressure
~~~~~~~~~~~~
Pressure of the ensemble, if required. Must be specified as a ``unyt_quantity``.

.. code-block:: python

  import unyt as u

  pressure=1*u.bar


n_molecules
~~~~~~~~~~~~~
List containing the number of molecules of each kind. Each elemnet of the list must be a positive integer.

thermostat
~~~~~~~~~~~~
Type of thermostat to be used. Must be a string. Default value is ``'NOSE'``.

traj_type
~~~~~~~~~~

Output trajectory format. Must be a string. Default is ``'PDB'``.

traj_freq
~~~~~~~~~~

Frequency of storing simulation trajectory. Must be positive integer.

seed
~~~~~
Seed for randomly initializing velocities. Must be a positive integer.

input_filename
~~~~~~~~~~~~~~
Name you want to give to the input file that is generated using cp2kmdpy. Default name will be depend on the project name.

output_filename
~~~~~~~~~~~~~~~~
Name you want to give to the output file. Default name will be depend on the project name.

initial_coordinate_filename
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In case structure generation from mBuild is not wanted/required, initial coordinates can be specified through this attribute.

use_atom_name_as_symbol
~~~~~~~~~~~~~~~~~~~~~~~~

If you want to make the atom name different from its symbol, set this flag to ``False``. In that case, you will have to set the ``particle.name`` and ``particle.element`` attributes of the particles in the ``mbuild`` molecule.

Running molecular dynamics
*********************************

0. Copy the setter.py file from ``cp2kmdpy/setter.py`` and change any settings if needed
1. Instantiate a ``cp2kmdpy.md.MD`` class and set all the attributes
2. Run the ``md_initialization`` method on that instance
3. Generate input files using ``setter.md_files`` function
4. Run using the ``run_md`` method from ``cp2kmdpy.runners``


Examples
*************

Dinitrogen *NpT* MD
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

  # ### Loading modules
  import numpy as np
  import unyt as u
  import mbuild as mb
  from cp2kmdpy.molecule_optimization import Molecule_optimization # for single molecule optimization
  from cp2kmdpy.md import MD # for running MD
  from cp2kmdpy import runners
  import setter


  # Defining the molecule we want to simulate


  class N2(mb.Compound):
      def __init__(self):
          super(N2, self).__init__()

          nitrogen1= mb.Particle(pos=[0.0, 0.0, 0.0], name='N')
          nitrogen2= mb.Particle(pos=[0.2, 0.0, 0.0], name='N')
          self.add([nitrogen2,nitrogen1])
          self.add_bond((nitrogen2,nitrogen1))

  #creating an instance of the MD class, also defining the parametrs of our md simulation
  molecule_list=[N2()]
  box=mb.box.Box(lengths=[1,1,1])
  q=MD(molecules=molecule_list,box=box,cutoff=200,functional='PBE',basis_set={'N':'DZVP-MOLOPT-GTH'},periodicity='XYZ',n_molecules=[2],traj_type='PDB',seed=1,project_name='N2_npt')


  #Setting temperature and ensemble, could have also been set while creating an instance of the MD class

  q.temperature=273.15*u.K
  q.ensemble='NPT_I'
  q.simulation_time=10*u.fs
  q.pressure=1*u.bar

  #Initializing q

  q.md_initialization()


  #generating input files
  setter.md_files(q)


  #running md
  runners.run_md(q.input_filename,q.output_filename )
