#!/usr/bin/env python3
# 
# Network spreading models using personalised (age/gender-matched) connectomes 
# from PPMI healthy controls for anatomically constrained tractography
#
# PPMI data wrangling: prerequisite for act_ppmi_1_preprocess.py
#
# Neil Oxtoby, February 2020

import os,subprocess,glob,csv
import pandas as pd
# import json

#* Name of CSV file containing image file paths and IDs: to be used as input to act_ppmi_1_preprocess.py
csvFile_Input = 'act_ppmi_filepaths.csv'
createCSV = True

#* File locations
top_folder = '~/ox_tools/anatomically_constrained_tractography'
data_folder = os.path.join(top_folder,'PPMI')
dicom_folder = os.path.join(data_folder,'raw_data_from_LONI')
nifti_folder = os.path.join(data_folder,'dicom_to_nifti')
sMRI_folder = os.path.normpath(os.path.join(data_folder,'act','sMRI'))
dMRI_folder = os.path.normpath(os.path.join(data_folder,'act','dMRI'))
json_folder = os.path.normpath(os.path.join(data_folder,'act','json'))
o = subprocess.run(['mkdir','-p',sMRI_folder])
o = subprocess.run(['mkdir','-p',dMRI_folder])
o = subprocess.run(['mkdir','-p',json_folder])

wd = os.path.join(top_folder,'scripts')
os.chdir(wd)
import act_utilities

#* Convert LONI DICOMs to NIFTI
dcm2niix_out = subprocess.run(['dcm2niix','-o',nifti_folder,'-f','%i_%t_%s_%p','-z','y',data_folder])

#* Move the nifti files in sMRI and dMRI folders (converted using dcm2niix from dicoms)
df_nifti = pd.DataFrame(columns=["Subject ID","Modality","Path"])

dmris = glob.glob(os.path.join(nifti_folder,"*DTI*"))
for dMRI in dmris:
    patno = os.path.basename(dMRI).split('_')[0]
    if '.json' in dMRI:
        o = subprocess.run(['mv',dMRI,os.path.normpath(json_folder)+"/"])
    else:
        o = subprocess.run(['mv',dMRI,os.path.normpath(dMRI_folder)+"/"])
        if ".nii" in dMRI:
            newpath = os.path.join(dMRI_folder,os.path.basename(dMRI))
            if len(df_nifti.index.values)==0:
                idx = 0
            else:
                idx = df_nifti.index.values[-1] + 1
            df_nifti.at[idx,'Subject ID'] = patno
            df_nifti.at[idx,'Modality']   = 'dMRI'
            df_nifti.at[idx,'Path']       = newpath
smris = glob.glob(os.path.join(nifti_folder,"*MPRAGE*"))
for sMRI in smris:
    patno = os.path.basename(sMRI).split('_')[0]
    if '.json' in sMRI:
        o = subprocess.run(['mv',sMRI,os.path.normpath(json_folder)+"/"])
    else:
        o = subprocess.run(['mv',sMRI,os.path.normpath(sMRI_folder)+"/"])
        newpath = os.path.join(sMRI_folder,os.path.basename(sMRI))
        if len(df_nifti.index.values)==0:
            idx = 0
        else:
            idx = df_nifti.index.values[-1] + 1
        df_nifti.at[idx,'Subject ID'] = patno
        df_nifti.at[idx,'Modality']   = 'sMRI'
        df_nifti.at[idx,'Path']       = newpath

if createCSV:
    df_nifti.sort_values(by=['Subject ID','Modality']).to_csv(os.path.join(top_folder,csvFile_Input),index=False)

#*** gzip the NIFTI files
# mri_ni = df_nifti["Path"].tolist()
# mri_gz = [n.replace('.nii','.nii.gz') for n in mri_ni]
# for j,k in zip(mri_ni,mri_gz):
#     if not os.path.exists(k):
#         oot = subprocess.call(['mrconvert',j,k])
#         oot2 = subprocess.call(['rm',j])
#     else:
#         print('File exists, moving along: {0}'.format(k))
# #* Update CSV
# if createCSV:
#     print('Writing CSV file: {0}'.format(csvFile_Input))
#     df_nifti['Path'].map(lambda x: x.replace('.nii','.nii.gz'))
#     df_nifti.sort_values(by=['Subject ID','Modality']).to_csv(os.path.join(top_folder,csvFile_Input),index=False)
# else:
#     print('Not creating CSV file {0}'.format(csvFile_Input))
