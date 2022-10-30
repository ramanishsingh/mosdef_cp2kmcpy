import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import mbuild as mb
import foyer
from foyer import Forcefield

import mbuild.formats.charmm_writer as mf_charmm

import mosdef_cp2k_writer


import parmed as pmd
from constrainmol import ConstrainedMolecule
import unyt as u

from mosdef_cp2kmcpy.mc import MC
import setter
warnings.filterwarnings("ignore", category=DeprecationWarning)
 

hydrogen_res_name = 'H2'
FF_file_hydrogen = 'h2_ff.xml'
hydrogen_FF = Forcefield(forcefield_files=FF_file_hydrogen)
hydrogen_mb = mb.load('[HH]', smiles=True)
for particle in hydrogen_mb.particles():
    particle.name = "HH"
    particle.element = "H"
    
    
hydrogen_mb.name = hydrogen_res_name

hydrogen_typed=hydrogen_FF.apply(hydrogen_mb)



constrain_mol = ConstrainedMolecule(hydrogen_typed)
constrain_mol.solve()
hydrogen_optimized_mb=mb.clone(hydrogen_mb)
hydrogen_optimized_mb.xyz=constrain_mol.xyz/10


water_res_name = 'H2O'
FF_file_water = 'water_ff.xml'
water_FF = Forcefield(forcefield_files=FF_file_water)
water_mb = mb.load('O', smiles=True)

for particle in water_mb.particles():
    
    if particle.name == "O":
        particle.element = "O"
    else:
        particle.element = "H"
        
water_mb.name = water_res_name

water_typed=water_FF.apply(water_mb)



constrain_mol = ConstrainedMolecule(water_typed)
constrain_mol.solve()
water_optimized_mb=mb.clone(water_mb)
water_optimized_mb.xyz=constrain_mol.xyz/10

#random box to be used in mosdef_charmm_writer
hydrogen_box_bias  = mb.fill_box(compound=[hydrogen_optimized_mb ],
                             n_compounds=[1] ,
                            box=[1.0, 1.0, 1.0] )

water_box_bias  = mb.fill_box(compound=[water_optimized_mb ],
                             n_compounds=[1] ,
                            box=[1.0, 1.0, 1.0] )

mixed_box  = mb.fill_box(compound=[hydrogen_optimized_mb, water_box_bias ],
                             n_compounds=[1, 1] ,
                            box=[1.0, 1.0, 1.0] )



FF_Dict = {hydrogen_optimized_mb.name:FF_file_hydrogen, water_optimized_mb.name:FF_file_water }

residues_List = [hydrogen_optimized_mb.name, water_optimized_mb.name ]




#hydrogen_box_bias.save("hydrogen_bias_coord.xyz",overwrite=True)
#nitrogen_box_bias.save("nitrogen_bias_coord.xyz",overwrite=True)

print('Running: GOMC FF file, and the psf and pdb files for the biasing potential file')
mf_charmm.charmm_psf_psb_FF(mixed_box ,
                            'overall_bias',
                            FF_filename ="overall_charmm_bias" ,
                            forcefield_selection = FF_Dict,
                            residues= residues_List ,
                            bead_to_atom_name_dict = {},
                            fix_residue = None,
                            reorder_res_in_pdb_psf = False
                            )


FF_Dict = {hydrogen_optimized_mb.name:FF_file_hydrogen }

residues_List = [hydrogen_optimized_mb.name ]




#hydrogen_box_bias.save("hydrogen_bias_coord.xyz",overwrite=True)
#nitrogen_box_bias.save("nitrogen_bias_coord.xyz",overwrite=True)

print('Running: GOMC FF file, and the psf and pdb files for the biasing potential file')
mf_charmm.charmm_psf_psb_FF(hydrogen_box_bias ,
                            'hydrogen_bias',
                            FF_filename ="hydrogen_charmm_bias" ,
                            forcefield_selection = FF_Dict,
                            residues= residues_List ,
                            bead_to_atom_name_dict = {},
                            fix_residue = None,
                            reorder_res_in_pdb_psf = False
                            )



FF_Dict = { water_optimized_mb.name:FF_file_water }

residues_List = [water_optimized_mb.name ]




#hydrogen_box_bias.save("hydrogen_bias_coord.xyz",overwrite=True)
#nitrogen_box_bias.save("nitrogen_bias_coord.xyz",overwrite=True)

print('Running: GOMC FF file, and the psf and pdb files for the biasing potential file')
mf_charmm.charmm_psf_psb_FF(water_box_bias ,
                            'water_bias',
                            FF_filename ="water_charmm_bias" ,
                            forcefield_selection = FF_Dict,
                            residues= residues_List ,
                            bead_to_atom_name_dict = {},
                            fix_residue = None,
                            reorder_res_in_pdb_psf = False
                            )



molecule_list=[hydrogen_optimized_mb, water_optimized_mb]
box=mb.box.Box(lengths=[1,1,1])
#Comment why pressure is needed

q=MC(molecule_list=molecule_list,n_box=2,n_molecules_each_box=[[2,2],[2,2]], box_list=[box,box],cutoff=200,functional='PBE',
     basis_set={'H':'DZVP-MOLOPT-GTH', 'O':'DZVP-MOLOPT-GTH','HH':'DZVP-MOLOPT-GTH'}, periodicity=['XYZ']*2,ensemble='GEMC_NPT',seed=1,project_name='water_NPT_GEMC',restart='FALSE',pressure=1*u.bar, use_atom_name_as_symbol = False)#, initial_coordinate_filename = ["box1_xyz.xyz", "box2_xyz.xyz"])


q.mc_initialization()


q.topology_filename=['hydrogen_bias.psf', 'water_bias.psf']

# move_probabilities=[pmavbmc,pmcltrans,pmhmc,pmswap,pmtraion,pmtrans,pmvolume]
#volume moves = PMVOLUME
# swap moves = PMSWAP - PMVOLUME
# AVBMC moves = PMAVBMC - PMSWAP - PMVOLUME
# an “inner” move = 1.0 - (PMAVBMC + PMSWAP + PMVOLUME)
#conformational changes = “inner” move percentage × PMTRAION
# molecular translation = “inner” move percentage × (PMTRANS - PMTRAION)
# molecular rotation = “inner” move percentage × (1.0 - PMTRANS - PMTRAION)

q.move_probabilities=[0,0.0,0,0.4,0.5,0.75,0.05]



# mol_probabilities=[[PMAVBMC_MOL,PMSWAP_MOL , PMTRAION_MOL, PMTRANS_MOL,PMROT_MOL],[]]
q.mol_probabilities=[[[1, 1],[1, 1],[1, 1],[1, 1],[1, 1]],[[1, 1],[1, 1],[1, 1],[1, 1],[1, 1]]]



#avbmc probabilities=[[AVBMC_ATOM,AVBMC_RMIN,AVBMC_RMAX,PBIAS],[]]
q.avbmc_probabilities=[[[1, 1],[1, 1],[1, 1],[1, 1]],[[1, 1],[1, 1],[1, 1],[1, 1]]]

q.charmm_potential_file='overall_charmm_bias.inp'

setter.mc_files(q)

