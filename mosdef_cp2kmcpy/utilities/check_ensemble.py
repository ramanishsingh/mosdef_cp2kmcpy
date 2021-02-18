def pressure_ensemble(val):
    """

    :param val: string, Name of the ensemble
    :return: boolean, returns True if pressure should be specified
    """
    if val=='NPE_F' or val=='NPE_I' or val=='NPT_F' or val=='NPT_I' or val=='NPT_GEMC':
        return True
    else:
        return False

def temperature_ensemble(val):
    """

    :param val: string, Name of the ensemble
    :return: boolean, returns True if temperature should be specified
    """
    if val=='MSST' or val=='MSST_DAMPED' or val=='NPT_F' or val=='NPT_I' or val=='NVT' or val=='NVT_ADIABATIC' or val=='NVT_GEMC' or val=='NPT_GEMC':
        return True
    else:
        return False
