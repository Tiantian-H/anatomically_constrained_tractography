#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=10:00:0
#$ -S /bin/bash
#$ -N csh_preproc
#$ -j y

hostname
date

source /share/apps/source_files/conda.source
conda activate lilypad_env #py2.7_mrtrix_env 
export PATH=${PATH}:/share/apps/mrtrix3/bin
export PATH=${PATH}:/share/apps/cmic/niftyreg_v1.5.43/bin
source /share/apps/source_files/cuda/cuda-10.1.source
export PATH=${PATH}:/share/apps/cmic/fsl-5.0.10/bin
FSLDIR=/share/apps/cmic/fsl-5.0.10
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH
export PATH=${PATH}:/share/apps/ants-2.2.0/bin
export PATH=${PATH}:/share/apps/freesurfer-6.0.0/bin
export FREESURFER_HOME=/share/apps/freesurfer-6.0.0
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=***your_directoy***/anatomically_constrained_tractography/ADNI/act

cd ***your_directoy***/anatomically_constrained_tractography/scripts

python act_adni_p1_preprocess.py
date
