# Deleting files that we do not need, files generated from a previous run
import os
import glob

extension_list=["*inp*","*out","*ener","*rest*","*Hess*","*REST*","*.xyz","*.pdb"]
for name in extension_list:
    
    filelist=glob.glob(name)
    for file in filelist:
        os.remove(file)


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
box=mb.box.Box(lengths=[1,1,1])


# First optimizng the structure of a single oxygen molecule

oxygen_optimization=Molecule_optimization(molecule=molecule,basis_set={'O':'DZVP-MOLOPT-GTH'},box=box,cutoff=300,functional='PBE',periodicity='NONE',n_iter=1)

#Initializing the optimization
oxygen_optimization.optimization_initialization()


#Generating optimization input files
setter.single_molecule_opt_files(oxygen_optimization)

#Running optimization
oxygen_optimization.run_optimization()


#Retrieving the structure of optimized molecule
optimized_oxygen=oxygen_optimization.return_optimized_molecule()


#checking the find bond length
np.linalg.norm(optimized_oxygen.xyz[1]-optimized_oxygen.xyz[0])


#creating an instance of the MD class, also defining the parametrs of our md simulation

q=MD(molecules=[optimized_oxygen],box=box,cutoff=200,functional='PBE',basis_set={'O':'DZVP-MOLOPT-GTH'},periodicity='XYZ',n_molecules=[2],traj_type='PDB',seed=1,project_name='twoO2molecules')


#Setting temperature and ensemble, could have also been set while creating an instance of the MD class

q.temperature=273.15*u.K
q.ensemble='NVT'
q.simulation_time=10*u.fs
#Initializing q
q.md_initialization()


#generating input files
setter.md_files(q)


#running md
runners.run_md(q.input_filename,q.output_filename )


