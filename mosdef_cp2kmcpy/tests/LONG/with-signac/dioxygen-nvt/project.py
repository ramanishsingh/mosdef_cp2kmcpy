import flow
from flow import FlowProject, directives
import warnings
import shutil
import os

warnings.filterwarnings("ignore", category=DeprecationWarning)


class Project(FlowProject):
    pass


@Project.label
def has_setter(job):
    return(job.isfile('setter.py'))



@Project.label
def has_md_files(job):
    "Verify that the md simulation has completed"

    try:
        return(job.isfile(job.doc.input_filename))
    except:
        return False

@Project.label
def md_completed(job):
    "Verify that the md simulation has completed"

    try:
        return(job.isfile(job.doc.output_filename))
    except:
        return False



@Project.operation
@Project.post(has_setter)
def copy_setter(job):
    with job:
        print(os.listdir())
        shutil.copyfile(Project().root_directory()+'/setter.py',job.workspace()+'/setter.py')



@Project.operation
@Project.post(has_md_files)
def md_files(job):
    with job:
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


        temperature = job.sp.T * u.K
        #Defining the molecule we want to simulate

        class O2(mb.Compound): # this class builds an oxygen molecule with a bond-length given in the oxygen2 x coor (nm)
            def __init__(self):
                super(O2, self).__init__()
        
                oxygen1= mb.Particle(pos=[0.0, 0.0, 0.0], name='O')
                oxygen2= mb.Particle(pos=[0.15, 0.0, 0.0], name='O')
                self.add([oxygen2,oxygen1])
                self.add_bond((oxygen2,oxygen1))


        molecule=O2();
        box=mb.box.Box(lengths=[1,1,1])
        q=MD(molecules=[molecule],box=box,cutoff=600,functional='PBE',basis_set={'O':'DZVP-MOLOPT-GTH'},periodicity='XYZ',n_molecules=[20],traj_type='PDB',seed=1,project_name='twentyO2molecules')
        q.temperature=temperature
        q.ensemble='NVT'
        q.simulation_time=1*u.ps
        #Initializing q
        q.md_initialization()


        #generating input files
        setter.md_files(q)
        job.doc.input_filename=q.input_filename
        job.doc.output_filename=q.output_filename

@Project.operation
@Project.pre(has_md_files)
@Project.post(md_completed)
@directives(np=4)
def run_md(job):
    from cp2kmdpy import runners
    with job:
        runners.run_md(job.doc.input_filename,job.doc.output_filename)

if __name__ == "__main__":
    Project().main()
