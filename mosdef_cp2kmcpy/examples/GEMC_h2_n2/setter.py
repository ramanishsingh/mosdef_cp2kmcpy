
import mbuild as mb
import subprocess
from mosdef_cp2k_writer.classes import SIM as sim
import ele
import numpy as np
from mosdef_cp2kmcpy.utilities import potential_setter
from mosdef_cp2kmcpy.utilities import set_basis_set
from mosdef_cp2kmcpy.utilities import check_ensemble




def mc_files(instance):
    molecule_list=instance.molecule_list
    n_box=instance.n_box
    n_molecules_each_box=instance.n_molecules_each_box
    box_list=instance.box_list
    initial_coordinate_filename=instance.initial_coordinate_filename
    use_atom_name_as_symbol=instance.use_atom_name_as_symbol
    cutoff = instance.cutoff
    scf_tolerance=instance.scf_tolerance
    basis_set_filename=instance.basis_set_filename
    basis_set=instance.basis_set
    potential_filename=instance.potential_filename
    functional=instance.functional
    periodicity=instance.periodicity
    n_steps=instance.n_steps
    n_ff_moves=instance.n_ff_moves
    nswapmoves=instance.nswapmoves
    ensemble=instance.ensemble
    project_name=instance.project_name
    temperature=instance.temperature
    pressure=instance.pressure
    traj_type=instance.traj_type
    traj_freq=instance.traj_freq
    seed=instance.seed
    input_filename=instance.input_filename
    output_filename=instance.output_filename
    move_probabilities=instance.move_probabilities
    mol_probabilities=instance.mol_probabilities
    avbmc_probabilities=instance.avbmc_probabilities
    run_type=instance.run_type
    restart=instance.restart
    restart_filename=instance.restart_filename
    topology_filename=instance.topology_filename
    charmm_potential_file=instance.charmm_potential_file



    if initial_coordinate_filename is None:
        initial_coordinate_filename=[]
        filled_boxes=[]

        for box_number in range(n_box):
            box= box_list[box_number]

            filled_box=mb.packing.fill_box(
                compound=molecule_list, n_compounds=n_molecules_each_box[box_number], box=box
            )
            box_file_name = project_name +"_box{}_".format(box_number+1)+"initial"+ ".xyz"
            initial_coordinate_filename.append(box_file_name)
            filled_box.save(box_file_name, overwrite="True")
            filled_boxes.append(filled_box)
            with open(box_file_name, "r") as fin:
                data = fin.read().splitlines(True)
            with open(box_file_name, "w") as fout:
                fout.writelines(data[2:])  # deleting first two lines
    else:
        for box_number in range(n_box):

            filled_box = mb.load(initial_coordinate_filename[box_number])
            box_file_name = project_name +"_box{}_".format(box_number+1)+"initial"+ ".xyz"
            filled_box.save(box_file_name, overwrite="True")
            with open(box_file_name, "r") as fin:
                data = fin.read().splitlines(True)
            with open(box_file_name, "w") as fout:
                fout.writelines(data[2:])  # deleting first two lines
    print("MC initial structure saved as {}".format(initial_coordinate_filename))
 
    atom_list = []
    mass_list = []
    symbol_list = []
    for i in range(len(molecule_list)):
        current_molecule = mb.clone(molecule_list[i])
        for particle in current_molecule.particles():
            atom_list.append(particle.name)
            if not (use_atom_name_as_symbol):
                symbol_list.append(particle.element.symbol)
            else:
                symbol_list.append(particle.name)
            if use_atom_name_as_symbol:
                mass_list.append(
                    ele.element_from_symbol("{}".format(particle.name)).mass
                )
            else:
                mass_list.append(
                    ele.element_from_symbol("{}".format(particle.element.symbol)).mass
                )

   
    unique_atom_list, indices = np.unique(atom_list, return_index=True)
    unique_atom_list = list(unique_atom_list)

    unique_symbol_list = [symbol_list[i] for i in indices]
    num_atoms = len(atom_list)
    num_unique_atoms = len(unique_atom_list)
    #writing box 1 file

    mySim = sim.SIM()
    # setting defaults

    mySim = sim.SIM()
    mySim.GLOBAL.PRINT_LEVEL="LOW"
    mySim.GLOBAL.RUN_TYPE = "MC"
    mySim.GLOBAL.PROJECT_NAME = project_name
    mySim.GLOBAL.PRINT_LEVEL = "LOW"
    mySim.GLOBAL.SEED = seed





    # FORCE EVAL SECTION
    mySim.FORCE_EVAL.METHOD = "QUICKSTEP"
    mySim.FORCE_EVAL.STRESS_TENSOR = "ANALYTICAL"

    mySim.FORCE_EVAL.DFT.BASIS_SET_FILE_NAME = basis_set_filename
    mySim.FORCE_EVAL.DFT.POTENTIAL_FILE_NAME = potential_filename
    mySim.FORCE_EVAL.DFT.CHARGE = 0
    mySim.FORCE_EVAL.DFT.MULTIPLICITY = 1
    mySim.FORCE_EVAL.DFT.MGRID.CUTOFF = cutoff
    mySim.FORCE_EVAL.DFT.MGRID.REL_CUTOFF = 60
    mySim.FORCE_EVAL.DFT.MGRID.NGRIDS = 4
    mySim.FORCE_EVAL.DFT.QS.METHOD = "GPW"
    mySim.FORCE_EVAL.DFT.QS.EPS_DEFAULT = 1e-8
    mySim.FORCE_EVAL.DFT.QS.EXTRAPOLATION = "ASPC"
    mySim.FORCE_EVAL.DFT.POISSON.PERIODIC = periodicity[0]
    mySim.FORCE_EVAL.DFT.PRINT.E_DENSITY_CUBE.SECTION_PARAMETERS = "OFF"
    mySim.FORCE_EVAL.DFT.SCF.SCF_GUESS = "ATOMIC"
    mySim.FORCE_EVAL.DFT.SCF.MAX_SCF = 20
    mySim.FORCE_EVAL.DFT.SCF.EPS_SCF = scf_tolerance
    mySim.FORCE_EVAL.DFT.SCF.SCF_GUESS="RESTART"

    mySim.FORCE_EVAL.DFT.SCF.OT.SECTION_PARAMETERS = ".TRUE."
    mySim.FORCE_EVAL.DFT.SCF.OT.PRECONDITIONER = "FULL_SINGLE_INVERSE"
    mySim.FORCE_EVAL.DFT.SCF.OT.MINIMIZER = "DIIS"
    mySim.FORCE_EVAL.DFT.SCF.OUTER_SCF.SECTION_PARAMETERS = ".TRUE."

    mySim.FORCE_EVAL.DFT.SCF.OUTER_SCF.MAX_SCF = 10
    mySim.FORCE_EVAL.DFT.SCF.OUTER_SCF.EPS_SCF = 1e-6
    mySim.FORCE_EVAL.DFT.SCF.PRINT.RESTART.SECTION_PARAMETERS = "OFF"
    mySim.FORCE_EVAL.DFT.SCF.PRINT.DM_RESTART_WRITE = ".TRUE."

    mySim.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.SECTION_PARAMETERS = functional
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.POTENTIAL_TYPE = "PAIR_POTENTIAL"
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.TYPE = "DFTD3"
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.PARAMETER_FILE_NAME = "dftd3.dat"
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.REFERENCE_FUNCTIONAL = functional
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.R_CUTOFF = 12
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.CALCULATE_C9_TERM =".TRUE."
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.REFERENCE_C9_TERM =".TRUE."
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.LONG_RANGE_CORRECTION =".TRUE."
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.EPS_CN =1e-6
    mySim.FORCE_EVAL.DFT.XC.XC_GRID.XC_DERIV="SPLINE2"
    mySim.FORCE_EVAL.DFT.XC.XC_GRID.XC_SMOOTH_RHO="NONE"


    mySim.FORCE_EVAL.SUBSYS.CELL.ABC = "{a} {b} {c}".format(
        a=10 * box_list[0].lengths[0], b=10 * box_list[0].lengths[1], c=10 * box_list[0].lengths[2]
    )
    mySim.FORCE_EVAL.SUBSYS.CELL.ALPHA_BETA_GAMMA = "{a} {b} {c}".format(
        a=box_list[0].angles[0], b=box_list[0].angles[1], c=box_list[0].angles[2]
    )
    mySim.FORCE_EVAL.SUBSYS.CELL.CELL_REF.ABC="{a} {b} {c}".format(
        a=10 * box_list[0].lengths[0], b=10 * box_list[0].lengths[1], c=10 * box_list[0].lengths[2]
    )

    mySim.FORCE_EVAL.SUBSYS.COORD.DEFAULT_KEYWORD = project_name +"_box1_"+"initial"+ ".xyz"
    mySim.FORCE_EVAL.SUBSYS.CELL.PERIODIC=periodicity[0]
    mySim.FORCE_EVAL.SUBSYS.CELL.SYMMETRY="CUBIC"
    mySim.FORCE_EVAL.SUBSYS.init_atoms(num_atoms)

    for i in range(num_unique_atoms):
        mySim.FORCE_EVAL.SUBSYS.KIND[i + 1].SECTION_PARAMETERS = unique_atom_list[i]
        if basis_set == [None]:

            mySim.FORCE_EVAL.SUBSYS.KIND[i + 1].BASIS_SET = set_basis_set.basis_set_setter(
                unique_atom_list[i]
            )
        else:
            mySim.FORCE_EVAL.SUBSYS.KIND[i + 1].BASIS_SET = basis_set[
                unique_atom_list[i]
            ]
        if not (use_atom_name_as_symbol):
            mySim.FORCE_EVAL.SUBSYS.KIND[i + 1].ELEMENT = unique_symbol_list[i]

        mySim.FORCE_EVAL.SUBSYS.KIND[i + 1].POTENTIAL = potential_setter.potential(
            unique_atom_list[i], functional
        )

    mySim.FORCE_EVAL.SUBSYS.TOPOLOGY.CONN_FILE_FORMAT="MOL_SET"
    num_unique_molecules=len(molecule_list)
    mySim.FORCE_EVAL.SUBSYS.TOPOLOGY.MOL_SET.init_molecules(num_unique_molecules)
    for i in range(num_unique_molecules):

        mySim.FORCE_EVAL.SUBSYS.TOPOLOGY.MOL_SET.MOLECULE[i+1].NMOL=n_molecules_each_box[0][i]
        mySim.FORCE_EVAL.SUBSYS.TOPOLOGY.MOL_SET.MOLECULE[i+1].CONN_FILE_NAME=topology_filename[i]

    # MOTION SECTION
    mySim.MOTION.MC.ENSEMBLE=ensemble
    mySim.MOTION.MC.TEMPERATURE =temperature
    mySim.MOTION.MC.PRESSURE=pressure
    mySim.MOTION.MC.IPRINT= 1
    mySim.MOTION.MC.LBIAS='TRUE'
    mySim.MOTION.MC.LSTOP='FALSE'
    mySim.MOTION.MC.NMOVES=8
    mySim.MOTION.MC.NSWAPMOVES =nswapmoves
    mySim.MOTION.MC.NSTEP= n_steps
    mySim.MOTION.MC.RESTART=restart
    mySim.MOTION.MC.BOX2_FILE_NAME=input_filename[1]
    mySim.MOTION.MC.RESTART_FILE_NAME="mc_restart_1"

    mySim.MOTION.MC.ETA= [0, 0]
    pmavbmc=move_probabilities[0]
    pmcltrans=move_probabilities[1]
    pmhmc=move_probabilities[2]
    pmswap=move_probabilities[3]
    pmtraion=move_probabilities[4]
    pmtrans=move_probabilities[5]
    pmvolume=move_probabilities[6]






    mySim.MOTION.MC.MOVE_PROBABILITIES.PMSWAP= pmswap
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMTRAION= pmtraion
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMCLTRANS=pmcltrans
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMTRANS= pmtrans
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMVOLUME= pmvolume
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMHMC= pmhmc
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMAVBMC= pmavbmc


    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMAVBMC_MOL=mol_probabilities[0][0]
    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMSWAP_MOL =mol_probabilities[0][1]
    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMTRAION_MOL= mol_probabilities[0][2]
    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMTRANS_MOL= mol_probabilities[0][3]
    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMROT_MOL =mol_probabilities[0][4]

    mySim.MOTION.MC.MOVE_PROBABILITIES.BOX_PROBABILITIES.PMVOL_BOX =1.0
    mySim.MOTION.MC.MOVE_PROBABILITIES.BOX_PROBABILITIES.PMHMC_BOX= 1.0
    if run_type=='equilibration':

        mySim.MOTION.MC.MOVE_UPDATES.IUPTRANS =20
        mySim.MOTION.MC.MOVE_UPDATES.IUPVOLUME= 20
    elif run_type=='production':
        mySim.MOTION.MC.MOVE_UPDATES.IUPTRANS =2000000
        mySim.MOTION.MC.MOVE_UPDATES.IUPVOLUME= 2000000
    else:
        print('run_type should either be equilibration or production')

    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMDIHEDRAL= [3.0]*len(molecule_list)
    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMANGLE= [3.0]*len(molecule_list)
    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMBOND =[0.074]*len(molecule_list)
    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMROT =[26.0]*len(molecule_list)
    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMTRANS=[ 0.38]*len(molecule_list)
    mySim.MOTION.MC.MAX_DISPLACEMENTS.BOX_DISPLACEMENTS.RMVOLUME =100.5*len(molecule_list)

    mySim.MOTION.MC.AVBMC.AVBMC_ATOM =avbmc_probabilities[0][0]
    mySim.MOTION.MC.AVBMC.AVBMC_RMIN =avbmc_probabilities[0][1]
    mySim.MOTION.MC.AVBMC.AVBMC_RMAX =avbmc_probabilities[0][2]
    mySim.MOTION.MC.AVBMC.PBIAS= avbmc_probabilities[0][3]

    mySim.write_changeLog(fn="md-changeLog.out")
    mySim.write_errorLog()
    mySim.write_inputFile(fn=input_filename[0])
    print("MC input file saved as {}".format(input_filename[0]))


#writing box 2 file

    mySim = sim.SIM()
    # setting defaults

    mySim = sim.SIM()
    mySim.GLOBAL.PRINT_LEVEL="LOW"
    mySim.GLOBAL.RUN_TYPE = "MC"
    mySim.GLOBAL.PROJECT_NAME = project_name
    mySim.GLOBAL.PRINT_LEVEL = "LOW"
    mySim.GLOBAL.SEED = seed





    # FORCE EVAL SECTION
    mySim.FORCE_EVAL.METHOD = "QUICKSTEP"
    mySim.FORCE_EVAL.STRESS_TENSOR = "ANALYTICAL"

    mySim.FORCE_EVAL.DFT.BASIS_SET_FILE_NAME = basis_set_filename
    mySim.FORCE_EVAL.DFT.POTENTIAL_FILE_NAME = potential_filename
    mySim.FORCE_EVAL.DFT.CHARGE = 0
    mySim.FORCE_EVAL.DFT.MULTIPLICITY = 1
    mySim.FORCE_EVAL.DFT.MGRID.CUTOFF = cutoff
    mySim.FORCE_EVAL.DFT.MGRID.REL_CUTOFF = 60
    mySim.FORCE_EVAL.DFT.MGRID.NGRIDS = 4
    mySim.FORCE_EVAL.DFT.QS.METHOD = "GPW"
    mySim.FORCE_EVAL.DFT.QS.EPS_DEFAULT = 1e-8
    mySim.FORCE_EVAL.DFT.QS.EXTRAPOLATION = "ASPC"
    mySim.FORCE_EVAL.DFT.POISSON.PERIODIC = periodicity[1]
    mySim.FORCE_EVAL.DFT.PRINT.E_DENSITY_CUBE.SECTION_PARAMETERS = "OFF"
    mySim.FORCE_EVAL.DFT.SCF.SCF_GUESS = "ATOMIC"
    mySim.FORCE_EVAL.DFT.SCF.MAX_SCF = 20
    mySim.FORCE_EVAL.DFT.SCF.EPS_SCF = scf_tolerance
    mySim.FORCE_EVAL.DFT.SCF.SCF_GUESS="RESTART"

    mySim.FORCE_EVAL.DFT.SCF.OT.SECTION_PARAMETERS = ".TRUE."
    mySim.FORCE_EVAL.DFT.SCF.OT.PRECONDITIONER = "FULL_SINGLE_INVERSE"
    mySim.FORCE_EVAL.DFT.SCF.OT.MINIMIZER = "DIIS"
    mySim.FORCE_EVAL.DFT.SCF.OUTER_SCF.SECTION_PARAMETERS = ".TRUE."

    mySim.FORCE_EVAL.DFT.SCF.OUTER_SCF.MAX_SCF = 10
    mySim.FORCE_EVAL.DFT.SCF.OUTER_SCF.EPS_SCF = 1e-6
    mySim.FORCE_EVAL.DFT.SCF.PRINT.RESTART.SECTION_PARAMETERS = "OFF"
    mySim.FORCE_EVAL.DFT.SCF.PRINT.DM_RESTART_WRITE = ".TRUE."

    mySim.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.SECTION_PARAMETERS = functional
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.POTENTIAL_TYPE = "PAIR_POTENTIAL"
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.TYPE = "DFTD3"
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.PARAMETER_FILE_NAME = "dftd3.dat"
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.REFERENCE_FUNCTIONAL = functional
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.R_CUTOFF = 12
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.CALCULATE_C9_TERM =".TRUE."
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.REFERENCE_C9_TERM =".TRUE."
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.LONG_RANGE_CORRECTION =".TRUE."
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.EPS_CN =1e-6
    mySim.FORCE_EVAL.DFT.XC.XC_GRID.XC_DERIV="SPLINE2"
    mySim.FORCE_EVAL.DFT.XC.XC_GRID.XC_SMOOTH_RHO="NONE"


    mySim.FORCE_EVAL.SUBSYS.CELL.ABC = "{a} {b} {c}".format(
        a=10 * box_list[1].lengths[0], b=10 * box_list[1].lengths[1], c=10 * box_list[1].lengths[2]
    )
    mySim.FORCE_EVAL.SUBSYS.CELL.ALPHA_BETA_GAMMA = "{a} {b} {c}".format(
        a=box_list[1].angles[0], b=box_list[1].angles[1], c=box_list[1].angles[2]
    )
    mySim.FORCE_EVAL.SUBSYS.CELL.CELL_REF.ABC="{a} {b} {c}".format(
        a=10 * box_list[1].lengths[0], b=10 * box_list[1].lengths[1], c=10 * box_list[1].lengths[2]
    )

    mySim.FORCE_EVAL.SUBSYS.COORD.DEFAULT_KEYWORD = project_name +"_box2_"+"initial"+ ".xyz"
    mySim.FORCE_EVAL.SUBSYS.CELL.PERIODIC=periodicity[1]
    mySim.FORCE_EVAL.SUBSYS.CELL.SYMMETRY="CUBIC"
    mySim.FORCE_EVAL.SUBSYS.init_atoms(num_atoms)

    for i in range(num_unique_atoms):
        mySim.FORCE_EVAL.SUBSYS.KIND[i + 1].SECTION_PARAMETERS = unique_atom_list[i]
        if basis_set == [None]:

            mySim.FORCE_EVAL.SUBSYS.KIND[i + 1].BASIS_SET = set_basis_set.basis_set_setter(
                unique_atom_list[i]
            )
        else:
            mySim.FORCE_EVAL.SUBSYS.KIND[i + 1].BASIS_SET = basis_set[
                unique_atom_list[i]
            ]
        if not (use_atom_name_as_symbol):
            mySim.FORCE_EVAL.SUBSYS.KIND[i + 1].ELEMENT = unique_symbol_list[i]

        mySim.FORCE_EVAL.SUBSYS.KIND[i + 1].POTENTIAL = potential_setter.potential(
            unique_atom_list[i], functional
        )

    mySim.FORCE_EVAL.SUBSYS.TOPOLOGY.CONN_FILE_FORMAT="MOL_SET"
    num_unique_molecules=len(molecule_list)
    mySim.FORCE_EVAL.SUBSYS.TOPOLOGY.MOL_SET.init_molecules(num_unique_molecules)
    for i in range(num_unique_molecules):

        mySim.FORCE_EVAL.SUBSYS.TOPOLOGY.MOL_SET.MOLECULE[i+1].NMOL=n_molecules_each_box[1][i]
        mySim.FORCE_EVAL.SUBSYS.TOPOLOGY.MOL_SET.MOLECULE[i+1].CONN_FILE_NAME=topology_filename[i]

    # MOTION SECTION
    mySim.MOTION.MC.ENSEMBLE=ensemble
    mySim.MOTION.MC.TEMPERATURE =temperature
    mySim.MOTION.MC.PRESSURE=pressure
    mySim.MOTION.MC.IPRINT= 1
    mySim.MOTION.MC.LBIAS='TRUE'
    mySim.MOTION.MC.LSTOP='FALSE'
    mySim.MOTION.MC.NMOVES=8
    mySim.MOTION.MC.NSWAPMOVES =nswapmoves
    mySim.MOTION.MC.NSTEP= n_steps
    mySim.MOTION.MC.RESTART=restart
    mySim.MOTION.MC.BOX2_FILE_NAME=input_filename[0]
    mySim.MOTION.MC.RESTART_FILE_NAME="mc_restart_2"

    mySim.MOTION.MC.ETA= [0.0, 0]
    pmavbmc=move_probabilities[0]
    pmcltrans=move_probabilities[1]
    pmhmc=move_probabilities[2]
    pmswap=move_probabilities[3]
    pmtraion=move_probabilities[4]
    pmtrans=move_probabilities[5]
    pmvolume=move_probabilities[6]






    mySim.MOTION.MC.MOVE_PROBABILITIES.PMSWAP= pmswap
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMTRAION= pmtraion
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMCLTRANS=pmcltrans
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMTRANS= pmtrans
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMVOLUME= pmvolume
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMHMC= pmhmc
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMAVBMC= pmavbmc


    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMAVBMC_MOL=mol_probabilities[1][0]
    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMSWAP_MOL =mol_probabilities[1][1]
    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMTRAION_MOL= mol_probabilities[1][2]
    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMTRANS_MOL= mol_probabilities[1][3]
    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMROT_MOL =mol_probabilities[1][4]

    mySim.MOTION.MC.MOVE_PROBABILITIES.BOX_PROBABILITIES.PMVOL_BOX =1.0
    mySim.MOTION.MC.MOVE_PROBABILITIES.BOX_PROBABILITIES.PMHMC_BOX= 1.0
    if run_type=='equilibration':

        mySim.MOTION.MC.MOVE_UPDATES.IUPTRANS =20
        mySim.MOTION.MC.MOVE_UPDATES.IUPVOLUME= 20
    elif run_type=='production':
        mySim.MOTION.MC.MOVE_UPDATES.IUPTRANS =2000000
        mySim.MOTION.MC.MOVE_UPDATES.IUPVOLUME= 2000000
    else:
        print('run_type should either be equilibration or production')

    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMDIHEDRAL= [3.0]*len(molecule_list)
    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMANGLE= [3.0]*len(molecule_list)
    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMBOND =[0.074]*len(molecule_list)
    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMROT =[26.0]*len(molecule_list)
    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMTRANS=[ 0.38]*len(molecule_list)
    mySim.MOTION.MC.MAX_DISPLACEMENTS.BOX_DISPLACEMENTS.RMVOLUME =100.5*len(molecule_list)

    mySim.MOTION.MC.AVBMC.AVBMC_ATOM =avbmc_probabilities[1][0]
    mySim.MOTION.MC.AVBMC.AVBMC_RMIN =avbmc_probabilities[1][1]
    mySim.MOTION.MC.AVBMC.AVBMC_RMAX =avbmc_probabilities[1][2]
    mySim.MOTION.MC.AVBMC.PBIAS= avbmc_probabilities[1][3]

    mySim.write_changeLog(fn="md-changeLog.out")
    mySim.write_errorLog()
    mySim.write_inputFile(fn=input_filename[1])
    print("MC input file saved as {}".format(input_filename[1]))


    #creating bias_template.inp file


    mySim = sim.SIM()
    mySim.GLOBAL.RUN_TYPE = "MC"
    mySim.GLOBAL.PROJECT_NAME="bias_template"
    mySim.GLOBAL.PRINT_LEVEL="LOW"

    mySim.MOTION.MC.ENSEMBLE="TRADITIONAL"
    mySim.MOTION.MC.TEMPERATURE =temperature
    mySim.MOTION.MC.PRESSURE=pressure
    mySim.MOTION.MC.IPRINT= 1
    mySim.MOTION.MC.LBIAS='TRUE'
    mySim.MOTION.MC.LSTOP='FALSE'
    mySim.MOTION.MC.NMOVES=8
    mySim.MOTION.MC.NSWAPMOVES =640
    mySim.MOTION.MC.NSTEP= 1
    mySim.MOTION.MC.RESTART='FALSE'

    mySim.MOTION.MC.ETA= 0.0
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMSWAP= 0.0
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMTRAION= 0.20
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMTRANS= 0.60
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMVOLUME= 0.02
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMHMC= 0.0
    mySim.MOTION.MC.MOVE_PROBABILITIES.PMAVBMC= 0.0
    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMAVBMC_MOL=[1]
    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMSWAP_MOL =[1.0]
    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMTRAION_MOL= [1.0]
    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMTRANS_MOL= [1.0]
    mySim.MOTION.MC.MOVE_PROBABILITIES.MOL_PROBABILITIES.PMROT_MOL =[1.0]
    mySim.MOTION.MC.MOVE_PROBABILITIES.BOX_PROBABILITIES.PMVOL_BOX =1.0
    mySim.MOTION.MC.MOVE_PROBABILITIES.BOX_PROBABILITIES.PMHMC_BOX= 1.0
    mySim.MOTION.MC.MOVE_UPDATES.IUPTRANS =100
    mySim.MOTION.MC.MOVE_UPDATES.IUPVOLUME= 100
    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMDIHEDRAL= [3.0]
    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMANGLE= [3.0]
    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMBOND =[0.074]
    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMROT =[26.0]
    mySim.MOTION.MC.MAX_DISPLACEMENTS.MOL_DISPLACEMENTS.RMTRANS=[ 0.38]
    mySim.MOTION.MC.MAX_DISPLACEMENTS.BOX_DISPLACEMENTS.RMVOLUME =100.5
    mySim.MOTION.MC.AVBMC.AVBMC_ATOM =[1]
    mySim.MOTION.MC.AVBMC.AVBMC_RMIN =[1.0]
    mySim.MOTION.MC.AVBMC.AVBMC_RMAX =[5.0]
    mySim.MOTION.MC.AVBMC.PBIAS= [0.5]

    mySim.FORCE_EVAL.METHOD="FIST"
    mySim.FORCE_EVAL.MM.FORCEFIELD.PARMTYPE="CHM"
    mySim.FORCE_EVAL.MM.FORCEFIELD.PARM_FILE_NAME=charmm_potential_file
    mySim.FORCE_EVAL.MM.FORCEFIELD.SPLINE.EMAX_SPLINE=10000
    mySim.FORCE_EVAL.MM.POISSON.EWALD.EWALD_TYPE="EWALD"
    mySim.FORCE_EVAL.MM.POISSON.EWALD.EWALD_ACCURACY=1e-6
    mySim.FORCE_EVAL.MM.POISSON.EWALD.GMAX=25

    mySim.FORCE_EVAL.SUBSYS.CELL.ABC="10 10 10"
    mySim.FORCE_EVAL.SUBSYS.CELL.PERIODIC=periodicity[0]
    mySim.FORCE_EVAL.SUBSYS.CELL.SYMMETRY="CUBIC"
    mySim.FORCE_EVAL.SUBSYS.COORD.DEFAULT_KEYWORD="bias_coord.xyz"
    mySim.FORCE_EVAL.SUBSYS.TOPOLOGY.CONN_FILE_FORMAT="MOL_SET"
    mySim.FORCE_EVAL.SUBSYS.TOPOLOGY.MOL_SET.init_molecules(num_unique_molecules)
    for i in range(num_unique_molecules):

        mySim.FORCE_EVAL.SUBSYS.TOPOLOGY.MOL_SET.MOLECULE[i+1].NMOL=n_molecules_each_box[0][i]
        mySim.FORCE_EVAL.SUBSYS.TOPOLOGY.MOL_SET.MOLECULE[i+1].CONN_FILE_NAME=topology_filename[i]

    mySim.write_changeLog(fn="changeLog.out")
    mySim.write_errorLog()
    mySim.write_inputFile(fn='bias_template.inp')
    
