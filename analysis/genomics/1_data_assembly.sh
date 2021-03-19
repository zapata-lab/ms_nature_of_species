#$ -cwd
#$ -m ae
#$ -l h_rt=6:00:00,h_data=128G,highmem


source /u/local/Modules/default/init/modules.sh
module load anaconda3
source $CONDA_DIR/etc/profile.d/conda.sh
conda activate ipyrad

# ipyrad commands
ipyrad -p ../../data/CladeI/params-CladeI_cladeSpecificDelimitation.txt -s 1234567
# when you are done with ipyrad
conda deactivate



