def info_molecule(molecule):
    """

    :param molecule: mBuild molecule
    :return:  names of elements in the molecule (list of strings), mass of the elements (list of floats)
    """
    num_atoms=len(molecule.atoms)
    name_atoms=[];
    mass_atoms=[];
    for i in range(num_atoms):
        name_atoms.append(molecule.atoms[i].element_name)
        mass_atoms.append(molecule.atoms[i].mass)

    return name_atoms,mass_atoms
