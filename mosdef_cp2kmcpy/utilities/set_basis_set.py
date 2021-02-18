def basis_set_setter(element_symbol):
    """

    :param element_symbol: string, element symbol
    :return: string, Basis set name
    """
    #HERE I should define all the elements and all the basis set
    if element_symbol=='H':
        return "TZV2PX-MOLOPT-GTH"
    elif element_symbol=='F':
        return "TZV2PX-MOLOPT-GTH"
    elif element_symbol=='Cl':
        return "TZV2PX-MOLOPT-GTH"
    elif element_symbol=='I':
        return "DZVP-MOLOPT-SR-GTH"
