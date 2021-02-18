from subprocess import call
def run_single_molecule_optimization(input_filename,output_filename):
    
    call("cp2k.popt -i {} -o {}".format(input_filename,output_filename),shell=True)
    print('Input file name given to runner is {}'.format(input_filename))
    print('Output file name  is {}'.format(output_filename))



def run_md(input_filename,output_filename):
    call("cp2k.popt -i {} -o {}".format(input_filename,output_filename),shell=True)
    print('Input file name given to runner is {}'.format(input_filename))
    print('Output file name  is {}'.format(output_filename))
