import os
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()



setuptools.setup(
    name="mosdef_cp2kmcpy", # Replace with your own username
    version="0.0.1",
    author="Ramanish Singh",
    author_email="singh891@umn.edu",
    description='A Python input/output interface to CP2K, developed as part of the MoSDeF simulation suite.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/ramanishsingh/mosdef_cp2kmcpy',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
