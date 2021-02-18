import signac
import numpy as np
import unyt as u


def init_project():

    # Initialize project
    project = signac.init_project("nvt1x1x1")

    # Define temperature
    temperature = 298.0 * u.K

    state_point = {"T": float(temperature.in_units(u.K).value)}

    job = project.open_job(state_point)
    job.init()


if __name__ == "__main__":
    init_project()
