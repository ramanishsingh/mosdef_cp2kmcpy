def is_cubic(box):
    """

    :param box: mBuild box
    :return: boolean, True if the box is cubic, False otherwise
    """
    angles=box.angles
    for angle in angles:
        if math.isclose(angle, 90.0):
            continue
        else:
            return False
            break
    lengths=box.lengths
    a=lengths[0]
    for length in lengths:
        if math.isclose(length, a):
            continue
        else:
            return False
            break
    return True
