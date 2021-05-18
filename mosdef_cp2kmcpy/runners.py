from subprocess import call
def run_single_molecule_optimization(input_filename,output_filename):
    
    call("cp2k.popt -i {} -o {}".format(input_filename,output_filename),shell=True)
    print('Input file name given to runner is {}'.format(input_filename))
    print('Output file name  is {}'.format(output_filename))



def run_mc(input_filename,output_filename):
    print('Input file name given to runner is {}'.format(input_filename))
    print('Output file name  is {}'.format(output_filename))
    call("cp2k.popt -i {} -i {} -o {}".format(input_filename[0],input_filename[1],output_filename),shell=True)
