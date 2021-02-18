.. # define a hard line break for HTML
.. |br| raw:: html

   <br />

CP2K installation
=================
CP2K is a Fortran-based software for atomistic simulations. More information `here <https://www.cp2k.org/>`_.
CP2K can be installed using the following methods (for detailed installation instructions, visit `this <https://www.cp2k.org/howto:compile>`_ page):

Installation using package manager
-----------------------------------
On Ubuntu (Linux) systems, CP2K can be easily installed using the following command:

>>> sudo apt install cp2k

After the installation, you can check your CP2K version using the following command:

>>> cp2k.popt --version

This installation method is fairly easy. However, CP2K will not run with full efficiency with this type of installation as the package manager may not fetch the latest version of the code.

Installation from source
-------------------------
CP2K can be installed from source using the following commands (assuming you have an intel compiler, you'd need to change the commands if using a gnu compiler):

>>> module load intel
>>> git clone https://github.com/cp2k/cp2k.git
>>> cd cp2k
>>> git submodule update --init --recursive

Now you will have to create an arch file according to the specifications of your system (.sopt/.popt file). In my case I have a .popt file named "Linux-x86-64-intel.popt" in the arch folder.

>>> make -j 4 ARCH=Linux-x86-64-intel VERSION=popt

You can add the CP2K executable path to the .bashrc file so that the program can be called easily.


>>> LOC_CP2K="$(pwd)"
>>> export PATH="${LOC_CP2K}:$PATH"

Or you can directly add the cp2k executable path to your .bashrc

>>> cd exe/Linux-x86-64-intel
>>> pwd

Copy the path and add the following line to the .bashrc file : |br|
 export PATH="$PATH:~/test-cp2k/cp2k/exe/Linux-x86-64-intel/"  |br|
 (replace the string "~/test-cp2k/cp2k/exe/Linux-x86-64-intel/" with your path)

>>> cp2k.popt --version



