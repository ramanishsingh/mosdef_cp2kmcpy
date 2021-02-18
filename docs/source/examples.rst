Examples
========

If you need help setting up any particular type of simulation with ``cp2kmdpy``, you can e-mail me at singh891@umn.edu.

Dioxygen *NVT* MD
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


  class O2(mb.Compound):
      def __init__(self):
          super(O2, self).__init__()

          oxygen1= mb.Particle(pos=[0.0, 0.0, 0.0], name='Od',element='O')
          oxygen2= mb.Particle(pos=[0.1208, 0.0, 0.0], name='Od',element='O')
          self.add([oxygen1,oxygen2])
          self.add_bond((oxygen1,oxygen2))

  #creating an instance of the MD class, also defining the parametrs of our md simulation
  molecule_list=[O2()]
  box=mb.box.Box(lengths=[1,1,1])
  q=MD(molecules=molecule_list,box=box,cutoff=200,functional='PBE',basis_set={'N':'DZVP-MOLOPT-GTH'},periodicity='XYZ',n_molecules=[2],traj_type='PDB',seed=1,project_name='O2_nvt',use_atom_name_as_symbol=False)


  #Setting temperature and ensemble, could have also been set while creating an instance of the MD class

  q.temperature=273.15*u.K
  q.ensemble='NVT'
  q.simulation_time=10*u.fs
  q.pressure=1*u.bar

  #Initializing q

  q.md_initialization()


  #generating input files
  setter.md_files(q)


  #running md
  runners.run_md(q.input_filename,q.output_filename )

