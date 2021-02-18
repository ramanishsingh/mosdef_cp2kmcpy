{% extends base_script %}
{% block project_header %}
#PBS -l walltime=00:50:00,mem=40gb

module purge
module load mkl
module load fftw
module load intel/cluster/2018
module load conda
source activate /home/siepmann/singh891/anaconda3/envs/mosdef36
date >> execution.log
{{ super() }}
{% endblock %}
