Molecule optimization
======================

This functionality can be used to first optimize the structure of the molecules in the simulation environment.
Optimization of a molecule structure is optional before the main molecular dynamics run, but is generally recommended.


Molecule_optimization class
****************************

.. autoclass:: cp2kmdpy.molecule_optimization.Molecule_optimization
    :members:

Instantiating ``Molecule_optimization`` first requires the specification of the molecule, functional, box , and the basis set to be used.
It can be done as follows:

.. code-block:: python

  import cp2kmdpy
  mol_opt = cp2kmdpy.molecule_optimization.Molecule_optimization(
      molecule=molecule,
      functional=functional,
      box=box,
      basis_set=basis_set
  )

molecule
~~~~~~~~~
The ``molecule`` is an instance of ``mBuild.Compund``. For example, if the molecule is an ethane molecule it can be created as follows:

.. code-block:: python

  import mbuild

  molecule = mbuild.load("CC", smiles=True)

functional
~~~~~~~~~~~~
``functional`` is the name of the DFT XC functional to be used for the molecule optimization.
``functional`` should be a string.

.. code-block:: python

  functional = 'PBE'


box
~~~~
In CP2K, a simulation box needs to be defined for placing the molecule while optimization.
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
During optimization, sometimes constraining the movement of some atoms is desired. This can be achieved by specifying ``fixed_atoms``.
For example, if atom number 1 to 80 needs to be fixed, it can be specified as follows

.. code-block:: python

  fixed_list='1..80'

periodicity
~~~~~~~~~~~
This attribute controls the periodicity of the box. Default value is ``'XYZ'``
Consider a system which is periodic in only X and Y direction. It can be specified as follows:

.. code-block:: python

  periodicity='XY'

n_iter
~~~~~~~
Number of iterations to conducted. Should be a positive integer.

Running molecule optimization
*********************************

0. Copy the setter.py file from ``cp2kmdpy/setter.py`` and change any settings if needed
1. Instantiate a ``cp2kmdpy.molecule_optimization.Molecule_optimization`` class and set all the attributes
2. Run the ``optimization_initialization`` method on that instance
3. Generate input files using ``setter.single_molecule_opt_files`` function
4. Run molecule optimization using the ``run_optimization()`` method

Example
***********************

Dioxygen structure optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

  # ### Loading modules
  import numpy as np
  import unyt as u
  import mbuild as mb
  from cp2kmdpy.molecule_optimization import Molecule_optimization # for single molecule optimization
  from cp2kmdpy.md import MD # for running MD
  from cp2kmdpy import runners
  import setter

  #Defining the molecule we want to simulate

  class O2(mb.Compound): # this class builds an oxygen molecule with a bond-length given in the oxygen2 x coor (nm)
      def __init__(self):
          super(O2, self).__init__()

          oxygen1= mb.Particle(pos=[0.0, 0.0, 0.0], name='O')
          oxygen2= mb.Particle(pos=[0.15, 0.0, 0.0], name='O')
          self.add([oxygen2,oxygen1])
          self.add_bond((oxygen2,oxygen1))

  molecule=O2();

  #Defining an mBuild box
  box=mb.box.Box(lengths=[2,2,2])

  #Creating an instance of the molecule optimization class
  oxygen_optimization=Molecule_optimization(molecule=molecule,basis_set={'O':'DZVP-MOLOPT-GTH'},box=box,cutoff=600,functional='PBE',periodicity='NONE',n_iter=100)

  #Initializing the optimization
  oxygen_optimization.optimization_initialization()


  #Generating optimization input files
  setter.single_molecule_opt_files(oxygen_optimization)

  #Running optimization
  oxygen_optimization.run_optimization()


  #Retrieving the structure of optimized molecule
  optimized_oxygen=oxygen_optimization.return_optimized_molecule()

