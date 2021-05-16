import mbuild as mb
import unyt as u
import ele



def remove_duplicate(x):
    return list(dict.fromkeys(x))


class biasing_potential():
    r"""

    Class that makes the biasing potential file (bias_template.inp) for MC simulations in CP2K

    :param molecules: Molecule(s) in the simulation box
    :type molecule: list, each element is an mBuild molecule
    """

    def __init__(self,charge_list=[], bond_list=[], bend_list=[], lennard_jones_list=[]
                 ):

        self.charge_list=charge_list
        self.bond_list=bond_list
        self.bend_list=bend_list
        self.lennard_jones_list=lennard_jones_list


