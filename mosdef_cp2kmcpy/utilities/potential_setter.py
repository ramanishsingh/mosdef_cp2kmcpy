def potential(element_symbol,functional):
    """

    :param element_symbol: string, symbol of the element
    :param functional: string, functional name
    :return: string, GTH potential name to be used in the simulation
    """
    return "GTH-"+functional
