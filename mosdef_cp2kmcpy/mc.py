import mbuild as mb
import unyt as u
import ele



def remove_duplicate(x):
    return list(dict.fromkeys(x))


class MC():
    r"""

    Base class for running Monte Carlo simulation.

    :param molecule_list: Molecules in the system
    :type molecule: list, each element is an mBuild molecule
    :param n_box: Number of boxes
    :type n_box: int, optional, default is 2
    :param n_molecules_each_box: number of molecules of each type in the boxes
    :type n_molecules_each_box: list, each element in the list is an ordered list containing the number of molecules of each type in the box corresponding to that element
    :param box_list: Simulation boxes
    :type box_list: list, each element is an mbuild box
    :param initial_coordinate_filename: Names of the files containing initial coordinates
    :type initial_coordinate_filename: list, optional, each element is a string, length must be equal to n_box
    :param use_atom_name_as_symbol:  If you want the atom name to be same as its symbol. Useful if you have want to use different basis sets for a single atom type in different environment. If you are setting it to false, make sure that the particles in your mbuild molecule have particle.element and particle.name attributes set (particle.element will be the symbol and particle.name will be the name used).
    :type use_atom_name_as_symbol: boolean, optional
    :param cutoff: Plane wave cutoff (Ry) for DFT calculation
    :type cutoff: float, optional
    :param scf_tolerance: Tolerance for each SCF cycle
    :type scf_tolerance: float, optional
    :param basis_set: Basis set for each atomic kind
    :type basis_set: dictionary with key being the atomic symbol and value being the basis set, both strings
    :param basis_set_filename: Filename for the basis set
    :type basis_set_filename: string, optional, defaults to BASIS_MOLOPT
    :param potential_filename: Filename for the pseudopotential to be used
    :type potential_filename: string, optional, defaults to GTH_POTENTIALS
    :param functional: DFT XC functional to be used
    :type functional: string
    :param periodicity: Periodicity of the box
    :type periodiicity: string, optional, defaults to 'XYZ'
    :param n_steps: Number of Monte Carlo cycles
    :type  n_steps: int
    :param n_ff_moves: Number of force-field-based moves between first-principles steps
    :type  n_ff_moves: int, optional, deault value is 8
    :param nswapmoves: Number insertions to try during each swap move
    :type  nswapmvoes: int, optional
    :param ensemble: Simulation ensemble
    :type ensemble: string
    :param project_name: Name of the project
    :type project_name: string, optional
    :param temperature: Simulation temperature
    :type temperature: unyt quantity
    :param pressure: Simulation pressure
    :type pressure: unyt quantity
    :param traj_type: Output trajectory format
    :type traj_type: string, optional
    :param seed: Random number seed for MD simulation
    :type seed: integer, optional
    :param input_filename: Name(s) that should be given to the input file(s)
    :type input_filename: list of string(s), optional, size of the list should equal the total number of boxes
    :param output_filename: Name(s) that should be given to the output file(s)
    :type output_filename: list of string(s), optional, size of the list should equal the total number of boxes
    :param move_probabilities: probabilities of different moves, in this order:[pmavbmc,pmcltrans,pmhmc,pmswap,pmtraion,pmtrans,pmvolume]
    :type move_probabilities: list of float
    :param mol_probabilities: move probabilities for each molecule in all the boes
    :type mol_probabilities: list in which each element contains the mol_probabilities for the molecules in the box corresponding to that element. The order is [[PMAVBMC_MOL,PMSWAP_MOL , PMTRAION_MOL, PMTRANS_MOL,PMROT_MOL],[similarly for box2]].Here PMAVBMC_MOL, PMSWAP_MOL, etc. will also be lists defining the move probabilities of each molecule type. Example, if you have 2 molecules and two boxes, and you want to give equal move probabilities to both molecules then mol_probabilities can be defined as mol_probabilities=[[[0.5,1],[0.5,1],[0.5,1],[0.5,1],[0.5,1]],[[0.5,1],[0.5,1],[0.5,1],[0.5,1],[0.5,1]]].
    :param avbmc_probabilities: AVBMC move probabilities for each molecule in all the boxes
    :type avbmc_probabilities: Defined in a similar manner as avbmc_probabilities.
    :param run_type: Type of run, equilibration or production
    :type run_type: string, optional, default is "equilibration"
    :param restart: If you want to restart a simulation
    :type restart: boolean, optional, default os False
    :param restart_filename: Name(s) of the restart files
    :type restart_filename: list of strings, length should equal the number of boxes
    :param topology_filename: Names of the topology files
    :type topology_filename: list of strings, length should equal the total number of unique molecules
    :param charmm_potential_file: Name of the charmm potential file to be used a biasing potential, if not available, bias_template.inp file must be prepared
    :param charmm_potential_file: string
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
