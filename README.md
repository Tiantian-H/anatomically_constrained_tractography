# Anatomically Constrained Tractography for Neurodegenerative Disease Progression Modelling

[Neil Oxtoby](https://github.com/noxtoby), UCL Centre for Medical Image Computing
February 2022

## Further documentation on the ADNI connectome

Modified by Tiantian He,May 2022

The pakcage is modified from Neil Oxtoby's pakcage for generating connetomes using PPMI dataset: https://github.com/noxtoby/anatomically_constrained_tractography. It is written in python 3.
The original PPMI pakcage is modified in order to generate connetomes using the ADNI dataset.  

This documentation is a supplementary version to Neil's original documentation. It has been implemented on UCL cluster.

### Preliminaries


Please make sure the dcm2niix package has been installed to your conda environment. If not, please run:
```bash
conda install -c conda-forge dcm2niix 
```
And also make sure your have the following pakcages installed for image processing:
MRtrix3, Niftyreg, FSL, ANTS, Freesurfer.
If you're running the code on the UCL cluster, you can directly use the pre-installed pakcages on the cluster using:

Please replace ```***your_directoy***``` with your own directly on the cluster. 

```bash
source /share/apps/source_files/conda.source
conda activate lilypad_env 
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
export SUBJECTS_DIR= ***your_directoy***/anatomically_constrained_tractography

cd ***your_directoy***/anatomically_constrained_tractography/scripts
```
### File structure

The original file structure is:
```bash
***your_directoy**/anatomically_constrained_tractography
|-- act_adni_filepaths.csv
|-- ADNI
|   `-- raw_data_from_LONI
|       `-- ADNI
|           `-- 127_S_2234
|               |-- DTI
|               |   `-- I207971
|               |       |-- ADNI_127_S_2234_MR_Axial_DTI__br_raw_20101209083912264_224_S96119_I207971.dcm
|               |       |-- ADNI_127_S_2234_MR_Axial_DTI__br_raw_20101209083912736_1886_S96119_I207971.dcm
|               |       |-- ADNI_127_S_2234_MR_Axial_DTI__br_raw_20101209083912799_327_S96119_I207971.dcm
|               |       |-- ADNI_127_S_2234_MR_Axial_DTI__br_raw_20101209083912852_1844_S96119_I207971.dcm
|               |       |-- ADNI_127_S_2234_MR_Axial_DTI__br_raw_20101209083912949_532_S96119_I207971.dcm
|               |       ... other similar .dcm files for DTI images
|               `-- sMRI
|                   `-- I207966
|                       |-- ADNI_127_S_2234_MR_Sag_IR-SPGR__br_raw_20101209083912900_76_S96115_I207966.dcm
|                       |-- ADNI_127_S_2234_MR_Sag_IR-SPGR__br_raw_20101209083913389_153_S96115_I207966.dcm
|                       |-- ADNI_127_S_2234_MR_Sag_IR-SPGR__br_raw_20101209083913434_147_S96115_I207966.dcm
|                       |-- ADNI_127_S_2234_MR_Sag_IR-SPGR__br_raw_20101209083913533_131_S96115_I207966.dcm
|                       |-- ... other similar .dcm files for sMRI images
```
The final file structure has been put at the end.

### Step 1

Please change the ```top_folder``` inside the .py file to your own directory.
```bash
cd ***your_directoy***/anatomically_constrained_tractography/scripts
python act_adni_0_prep.py
```
### Step 2
Run act_adni_p1_preprocess.py file on the cluster by submitting the following task to the cluster. Please make sure that the directories in the following .sh file has been set to your own.
```bash
qsub run_preproc_p1.sh
```
### Step 3 
Run Freesurfer segmentation
See the file called newFS_S2234_2010_12_08_207966.sh in sMRI_reg an example sh file for your reference. Please make sure the directory inside the sh file is altered.

### Step 4 
Run act_adni_p2_preprocess.py file on the cluster by submitting the following task to the cluster. Please make sure that the directories in the following .sh file has been set to your own.
```bash
qsub run_preproc_p2.sh
```

### Step 5

Run the following .py file to generate another .sh file.

```bash
python adni_tckgen_cluster.py
```
### Step 6
Create a file (eg. called fivettgen_freesurfer) under ```***your_directoy***/anatomically_constrained_tractography/ADNI/act```. Move 5tt sMRI images and the preprocessed DTI images to this file. 

### Step 7
```
cd ***your_directoy***/anatomically_constrained_tractography/scripts
python adni_tckgen_cluster.py fivettgen_freesurfer freesurfer
```
After running the above command, a new command will be generated automatically as an output, such as:
```bash
qsub ***your_directoy***/anatomically_constrained_tractography/ADNI/act/fivettgen_freesurfer/fivettgenfreesurfer_filelist_20220518.txt.sh
```

The final file structure is:
```bash
***your_directoy**/anatomically_constrained_tractography
|-- act_adni_filepaths.csv
|-- ADNI
|   |-- act
|   |   |-- act_adni.csv
|   |   |-- dMRI
|   |   |   |-- 127_S_2234_20101208145414_7_ADNI2_GE_3T_20.0m4_8cha-DTI.bval
|   |   |   |-- 127_S_2234_20101208145414_7_ADNI2_GE_3T_20.0m4_8cha-DTI.bvec
|   |   |   |-- 127_S_2234_20101208145414_7_ADNI2_GE_3T_20.0m4_8cha-DTI.mif
|   |   |   `-- 127_S_2234_20101208145414_7_ADNI2_GE_3T_20.0m4_8cha-DTI.nii.gz
|   |   |-- dMRI_biascorrect
|   |   |   |-- 127_S_2234_20101208145414_7_ADNI2_GE_3T_20.0m4_8cha-DTI_preproc_biascorrect.mif
|   |   |   `-- 127_S_2234_20101208145414_7_ADNI2_GE_3T_20.0m4_8cha-DTI_preproc_biascorrect.nii.gz
|   |   |-- dMRI_denoised
|   |   |   `-- 127_S_2234_20101208145414_7_ADNI2_GE_3T_20.0m4_8cha-DTI_denoised.mif
|   |   |-- dMRI_fods
|   |   |   |-- 2234_wm_fods_tournier.mif
|   |   |   |-- 2234_wm_fods_tournier_normed.mif
|   |   |   `-- 2234_wm_response_tournier.txt
|   |   |-- dMRI_mask
|   |   |   `-- 127_S_2234_20101208145414_7_ADNI2_GE_3T_20.0m4_8cha-DTI_preproc_biascorrect_mask.mif
|   |   |-- dMRI_preproc
|   |   |   `-- 127_S_2234_20101208145414_7_ADNI2_GE_3T_20.0m4_8cha-DTI_preproc.mif
|   |   |-- fivettgen_freesurfer
|   |   |   |-- 127_S_2234_20101208145414-DTImask.mif
|   |   |   |-- 127_S_2234_20101208145414-DTI_preproc_biascorrect.mif
|   |   |   |-- 127_S_2234_20101208145414-RF.txt
|   |   |   |-- 127_S_2234_20101208145414-sMRI_reg_aladin_5tt_freesurfer.mif
|   |   |   |-- 127_S_2234_20101208145414-tckgen8M_freesurfer.tck
|   |   |   |-- 127_S_2234_20101208145414-tckgen8M_freesurfer_tckedit4M.tck
|   |   |   |-- 127_S_2234_20101208145414-WMFODs_tournier.mif
|   |   |   |-- fivettgenfreesurfer.e4663621.1
|   |   |   |-- fivettgenfreesurfer_filelist_20220518_done.txt
|   |   |   |-- fivettgenfreesurfer_filelist_20220518.txt
|   |   |   |-- fivettgenfreesurfer_filelist_20220518.txt.sh
|   |   |   `-- fivettgenfreesurfer.o4663621.1
|   |   |-- json
|   |   |   `-- 127_S_2234_20101208145414_7_ADNI2_GE_3T_20.0m4_8cha-DTI.json
|   |   |-- sMRI
|   |   |   `-- 127_S_2234_20101208145414_2_ADNI2_GE_3T_20.0m4_8cha-sMRI.nii.gz
|   |   |-- sMRI_5tt
|   |   |   `-- 127_S_2234_20101208145414_2_ADNI2_GE_3T_20.0m4_8cha-sMRI_reg_aladin_5tt_freesurfer.mif
|   |   |-- sMRI_freesurfer
|   |   |   `-- ADNI_S_2234_2010_12_08
|   |   |       `-- mri
|   |   |           |-- aparc.a2009s+aseg.mgz
|   |   |           |-- aparc+aseg.mgz
|   |   |           |-- aparc+aseg.nii.gz
|   |   |           |-- aparc.DKTatlas+aseg.mgz
|   |   |           |-- aseg.auto.mgz
|   |   |           ... other outputs of Freesurfer...
|   |   |-- sMRI_nodes
|   |   |   `-- ADNI_S_2234_2010_12_08_nodes_freesurfer.mif
|   |   `-- sMRI_reg
|   |       |-- 127_S_2234_20101208145414_2_ADNI2_GE_3T_20.0m4_8cha-sMRI_reg_aladin.nii.gz
|   |       |-- 127_S_2234_20101208145414_2_ADNI2_GE_3T_20.0m4_8cha-sMRI_reg_aladin_transform.txt
|   |       `-- newFS_S2234_2010_12_08_207966.sh
|   |-- dicom_to_nifti
|   |   `-- 127_S_2234_20101208145414_2_ADNI2_GE_3T_20.0m4_8cha-Sag_IR-SPGR.json
|   `-- raw_data_from_LONI
|       `-- ADNI
|           `-- 127_S_2234
|               |-- DTI
|               |   `-- I207971
|               |       |-- ADNI_127_S_2234_MR_Axial_DTI__br_raw_20101209083912264_224_S96119_I207971.dcm
|               |       |-- ADNI_127_S_2234_MR_Axial_DTI__br_raw_20101209083912736_1886_S96119_I207971.dcm
|               |       |-- ADNI_127_S_2234_MR_Axial_DTI__br_raw_20101209083912799_327_S96119_I207971.dcm
|               |       |-- ADNI_127_S_2234_MR_Axial_DTI__br_raw_20101209083912852_1844_S96119_I207971.dcm
|               |       |-- ADNI_127_S_2234_MR_Axial_DTI__br_raw_20101209083912949_532_S96119_I207971.dcm
|               |       ... other similar .dcm files for DTI images
|               `-- sMRI
|                   `-- I207966
|                       |-- ADNI_127_S_2234_MR_Sag_IR-SPGR__br_raw_20101209083912900_76_S96115_I207966.dcm
|                       |-- ADNI_127_S_2234_MR_Sag_IR-SPGR__br_raw_20101209083913389_153_S96115_I207966.dcm
|                       |-- ADNI_127_S_2234_MR_Sag_IR-SPGR__br_raw_20101209083913434_147_S96115_I207966.dcm
|                       |-- ADNI_127_S_2234_MR_Sag_IR-SPGR__br_raw_20101209083913533_131_S96115_I207966.dcm
|                       |-- ... other similar .dcm files for sMRI images
```
### History
This code was used to generate healthy elderly connectomes for network spreading models of pathology. 

It was develop using [PPMI](https://ppmi-info.org) data and the [MRTrix3 ACT pipeline](https://mrtrix.readthedocs.io/en/dev/quantitative_structural_connectivity/act.html) (see also [this ISMRM tutorial](https://mrtrix.readthedocs.io/en/dev/quantitative_structural_connectivity/ismrm_hcp_tutorial.html) on [HCP](http://www.humanconnectomeproject.org/) data).


### Contents

- `PPMI` folder: for data, including outputs
- `scripts` folder: the code
  - `act_ppmi_0_prep.py`: data wrangling
  - `act_ppmi_1_preprocess.py`: preprocessing steps to run locally
  - `ppmi_tckgen_cluster.py`: tractography to run on the UCL CMIC/CS cluster
  - `act_utilities.py`: utilities

PPMI imaging data was downloaded as raw ("Original") DICOM files. Will require tweaking for other datasets, such as ADNI.
