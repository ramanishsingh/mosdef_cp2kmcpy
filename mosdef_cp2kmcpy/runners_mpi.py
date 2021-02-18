from subprocess import Popen, PIPE
def run_single_molecule_optimization(input_filename,output_filename,np):
    
    print('Input file name given to runner is {}'.format(input_filename))
    print('Output file name  is {}'.format(output_filename))

    process=Popen("mpirun -n {} cp2k.popt -i {} -o {}".format(np,input_filename,output_filename),shell=True, universal_newlines=True,stdin=PIPE, stdout=PIPE, stderr=PIPE )
    output, error = process.communicate();
    print (output);

def run_md(input_filename,output_filename,np):
    print('Input file name given to runner is {}'.format(input_filename))
    print('Output file name  is {}'.format(output_filename))


    #process=Popen("mpirun -n {} ~/test-cp2k/cp2k/exe/Linux-x86-64-intel/cp2k.popt -i {} -o {}".format(np,input_filename,output_filename),shell=True, universal_newlines=True,stdin=PIPE, stdout=PIPE, stderr=PIPE )
    process=Popen("mpirun -n {} cp2k.popt -i {} -o {}".format(np, input_filename,output_filename),shell=True, universal_newlines=True,stdin=PIPE, stdout=PIPE, stderr=PIPE )
    output, error = process.communicate();
    print (output);
    return output,error
