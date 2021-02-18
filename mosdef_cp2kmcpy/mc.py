import mbuild as mb
import unyt as u
import ele



def remove_duplicate(x):
    return list(dict.fromkeys(x))


class MC():
    r"""

    Base class for running Monte Carlo simulation.

    :param molecules: Molecule(s) in the simulation box
    :type molecule: list, each element is an mBuild molecule
    """

    def __init__(self,molecule_list=None, n_box=2,n_molecules_each_box=None,box_list=None,
                 initial_coordinate_filename=None,use_atom_name_as_symbol=True,cutoff=None,
                 scf_tolerance=None, basis_set=None,basis_set_filename=None, potential_filename=None,
                 functional=None,
                 periodicity=None,n_steps=None,n_ff_moves=None, nswapmoves=None,ensemble=None,project_name=None,
                 temperature=None,pressure=None,traj_type=None,
                 traj_freq=None,seed=None,input_filename=None,output_filename=None,move_probabilities=None,
                 mol_probabilities=None,avbmc_probabilities=None, run_type='equilibration',restart=False,restart_filename=None,
                 topology_filename=None,charmm_potential_file=None,
                 ):

        self.molecule_list=molecule_list
        self.n_box=n_box
        self.n_molecules_each_box=n_molecules_each_box
        self.box_list=box_list
        self.initial_coordinate_filename=initial_coordinate_filename
        self.use_atom_name_as_symbol=use_atom_name_as_symbol
        self.cutoff=cutoff
        self.scf_tolerance=scf_tolerance
        self.basis_set_filename=basis_set_filename
        self.basis_set=basis_set
        self.potential_filename=potential_filename
        self.functional=functional
        self.periodicity=periodicity
        self.n_steps=n_steps
        self.n_ff_moves=n_ff_moves
        self.nswapmoves=nswapmoves
        self.ensemble=ensemble
        self.project_name=project_name
        self.temperature=temperature
        self.pressure=pressure
        self.traj_type=traj_type
        self.traj_freq=traj_freq
        self.seed=seed
        self.input_filename=input_filename
        self.output_filename=output_filename
        self.move_probabilities=move_probabilities
        self.mol_probabilities=mol_probabilities
        self.avbmc_probabilities=avbmc_probabilities
        self.run_type=run_type
        self.restart=restart
        self.restart_filename=restart_filename
        self.topology_filename=topology_filename
        self.charmm_potential_file=charmm_potential_file


    def mc_initialization(self):

        molecule_list=self.molecule_list
        n_box=self.n_box
        n_molecules_each_box=self.n_molecules_each_box
        box_list=self.box_list
        initial_coordinate_filename=self.initial_coordinate_filename
        use_atom_name_as_symbol=self.use_atom_name_as_symbol
        cutoff=self.cutoff
        scf_tolerance=self.scf_tolerance
        basis_set_filename=self.basis_set_filename
        basis_set=self.basis_set
        potential_filename=self.potential_filename
        functional=self.functional
        periodicity=self.periodicity
        n_steps=self.n_steps
        n_ff_moves=self.n_ff_moves
        nswapmoves=self.nswapmoves
        ensemble=self.ensemble
        project_name=self.project_name
        temperature=self.temperature
        pressure=self.pressure
        traj_type=self.traj_type
        traj_freq=self.traj_freq
        seed=self.seed
        input_filename=self.input_filename
        output_filename=self.output_filename
        move_probabilities=self.move_probabilities
        mol_probabilities=self.mol_probabilities
        avbmc_probabilities=self.avbmc_probabilities
        run_type=self.run_type
        restart=self.restart
        restart_filename=self.restart_filename
        topology_filename=self.topology_filename
        charmm_potential_file=self.charmm_potential_file


        if project_name==None:
            self.project_name='MC_sample_project'
            print('project_name not specified, set as MC_sample_project')
        if cutoff==None:
            self.cutoff=600;
            print('cutoff not specified, set as 600')
        if scf_tolerance==None:
            self.scf_tolerance=1e-6
            print('scf_tolerance not specified, set as 1e-6')
        if basis_set_filename==None:
            self.basis_set_filename='BASIS_MOLOPT'
            print('basis_set_filename not defined, set as BASIS_MOLOPT')
        if potential_filename==None:
            self.potential_filename='GTH_POTENTIALS'
            print('potential_filename not specified, set as GTH_POTENTIALS')
        if periodicity==None:
            self.periodicity=['XYZ']*n_box
            print('periodicity not specified, set as XYZ for {} boxes'.format(n_box))
        if n_steps==None:
            self.n_steps=1000
            print('n_steps not specified, set as 1000')
        if n_ff_moves==None:
            self.n_ff_moves=8
            print('n_ff_moves not specified, set as 8')

        if nswapmoves==None:
            self.nswapmoves=640
            print('nswapmoves not specified, set as 640, will be ignore if n_box<2')

        if ensemble==None:
            self.ensemble='GEMC_NVT'
            print('ensemble not specified, set as GEMC_NVT')

        if traj_type==None:
            self.traj_type='XYZ'
            print('output trajectory format set as XYZ')
        if traj_freq==None:
            self.traj_freq=10
        if seed == None:
            self.seed=0


        if self.input_filename==None:
            self.input_filename=[self.ensemble+'_box1.inp',self.ensemble+'_box2.inp']
            print('input_filename not specified, set as {}'.format(self.input_filename))

        if self.output_filename==None:
            self.output_filename=self.project_name+'_mc_output.out'
            print('output_filename not specified, set as {}'.format(self.output_filename))
        


        if self.temperature is not None:

            self.temperature=(temperature.to('K')).value
        else:
            self.temperature=298
            print('temperature not defined, set as 298 K')
        if self.pressure is not None:

            self.pressure=(pressure.to('bar')).value
        else:
            self.pressure=1
            print('pressure not defined, set as 1 bar. If ensemble=="NVT_GEMC" then pressure value will not be used')

        print('You can change default settings in setter.mc_files')
