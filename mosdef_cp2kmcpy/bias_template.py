import mbuild as mb
import unyt as u
import ele



def remove_duplicate(x):
    return list(dict.fromkeys(x))


class biasing_potential():
    r"""

    Class that makes the biasing potential file (bias_template.inp) for MC simulations in CP2K. Please refer to CHARGE, BOND, BEND, and NONBONDED sections here:https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD.html.

    :param charge_list: list of charge dictionary for each atom, only norzero charges are required
    :type charge_list: list of dict, keywords are "ATOM" (string) and "CHARGE" (float)
    :param bond_list: list of bond dictionary for each bond
    :type bond_list: list of dict, keywords are "ATOMS" (string), "CS" (float), "K" (float), "KIND" (string), and "R0" (float)
    :param bend_list: list of angle bending dictionary for each angle
    :type bend_list: list of dict, keywords are "ATOMS" (string), "CB" (float), "K" (float), "KBS12" (float), "KBS32" (float), "KIND" (string), "KSS" (float), "LEGENDRE" (string), "R012" (float), "R032" (float), and "THETA0" (float)
    :param lennard_jones_list: list of dict, lennard-jones parameters for each atom
    :type lennard_jones_list: list of dict, keywords are "ATOMS" (string), "EPSILON" (float), "RCUT" (float), "RMAX" (float), "RMIN" (float), and "SIGMA" (float)
    """

    def __init__(self,charge_list=[], bond_list=[], bend_list=[], lennard_jones_list=[]
                 ):

        self.charge_list=charge_list
        self.bond_list=bond_list
        self.bend_list=bend_list
        self.lennard_jones_list=lennard_jones_list


