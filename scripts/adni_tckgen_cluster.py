#!/usr/bin/env python
# encoding: utf-8
"""
ppmi_tckgen_cluster.py

Modified from adni_tckgen_cluster.py

March 2020
Copyright (c) 2017-2020 Neil Oxtoby.  All rights reserved.
"""

help_message = '''
1. Find a list of Structural-Connectome-related MRtrix image files in a given folder
    - *5TT*.mif 
    - [optionally] *WMFOD*.mif
2. See if there exists a corresponding tracks file (from tckgen or tcksift)
3. Resubmit as an array job to do the following:
    - [if no FOD file exists] generate response functions and WM FODs from DWI: dwi2response => dwi2fod
    - tckgen/tcksift

The list of images to be processed is written in a text file 
- one file name per line - to be used in a Sun Grid Engine array job that 
performs the tckgen and tsksift steps of the structural connectome pipeline.

    Usage:
	python ppmi_tckgen_cluster.py folderContaining5TTandWMFODs

'''
import fnmatch
import os
import shutil
import glob
import datetime
import sys

def traverseFolders(topFolder='/SAN/medic/Net_Mod_MS/HCPData/ADNIGIF3_SC', wildcard='*_5TT_gif3.mif', depth=1):
    """Traverses the file hierarchy within topFolder looking for files named *_5TT_gif3.mif 
    (or other supplied wildcard) and returns a list containing the paths, including filename, 
    for each file."""
    if depth==1:
        print('traverseFolders: depth==1 => globbing')
        imageFiles = glob.glob('{0}/{1}'.format(topFolder,wildcard))
        return imageFiles
    elif depth==2:
        print('traverseFolders: depth==2 => globbing')
        imageFiles = glob.glob('{0}/*/{1}'.format(topFolder,wildcard))
        return imageFiles
    else:
        print('traverseFolders: depth!=1 => os.walk')
        matches = []
        for root, dirnames, filenames in os.walk(topFolder):
            for filename in fnmatch.filter(filenames, wildcard):
                matches.append(os.path.join(root, filename))
        return matches
    
# note1 : manally move the dMRI.mif into the target file
# note2 : add mrtrix source through conda
# note3 : PREFIXFOD=${PREFIX5TT:0:26} change the index from 5:10 into 0:26 due to different names in ADNI from PPMI
# note4 : change preprocessed_biascorrected into _preproc_biascorrect because I used the for_each part instead
# note5 : use 'dwi2response tournier' instead
# note6 : modify dwi2fod
# note7 : to do: add mtnormalise for wmfod
# note8 : delete mrresize: command not exist in v3.0.3
# note9 : tckgen: [ERROR] unknown option "-number", change number to select
# fivett=${PREFIX5TT}_5TT_${algorithm}.mif to ${PREFIX5TT}
# tck=${PREFIXFOD}  PREFIX5TT to PREFIXFOD, same as tck_reduced
def qsubSC(workingDir, textFileListOfFilenames, jobName, algorithm):
    """For automating qsub submission of an array job to the CMIC cluster at UCL, 
    for processing a list of pre-processed medical images provided in a text file."""
    
    h_vmem = '7.4G'
    tmem = '7.4G'
    
    with open(textFileListOfFilenames) as f:
        nImages = sum(1 for _ in f)
    print('Found %i lines of text\n' % nImages)
    qsubText = """#$ -S /bin/bash
#$ -l hostname=burns*             # required: MRtrix build node
#$ -l h_rt=5:59:00               # Request wall time
#$ -l h_vmem=%s,tmem=%s          # Request memory
#$ -t 1-%s                        # Array job
#$ -N %s                   # $JOB_NAME
#$ -wd %s     # Working directory
# $ -l tscratch=15G
# $ -V # export environment variables
# $ -cwd # execute job from the current working directory
# $ -j y # merge stdout with stderr
# $ -R y # reservation y/n

echo "****** Job started ******"
date
START=$(date +%s)



tckgen_num=6M
tckedit_num=3M
tcksift_num=1M
algorithm=%s
rf_algo=tournier

source /share/apps/source_files/conda.source
conda activate lilypad_env 
export PATH=${PATH}:/share/apps/mrtrix3/bin

DATA_PATH=%s
FILELIST=%s
FILEPATH=$(awk "NR==$SGE_TASK_ID" $FILELIST)
PREFIX5TT=`basename ${FILEPATH} _5TT_${algorithm}.mif`
PREFIXFOD=${PREFIX5TT:0:26}
fivett=${PREFIX5TT}
wmfods=${PREFIXFOD}WMFODs_${rf_algo}.mif
wmfods_sft=${PREFIXFOD}WMFODs_${rf_algo}_downsampled.mif
tck=${PREFIXFOD}tckgen${tckgen_num}_${algorithm}.tck
tck_reduced=${PREFIXFOD}tckgen${tckgen_num}_${algorithm}_tckedit${tckedit_num}.tck
sft=`basename ${tck_reduced} .tck`downsampled_tcksift${tcksift_num}.tck

dwi=${PREFIXFOD}DTI_preproc_biascorrect.mif
dwimask=${PREFIXFOD}DTImask.mif
rf=${PREFIXFOD}RF.txt
rfvoxels=${PREFIXFOD}RFvoxels${algorithm}.mif

###### Generate DTI mask
if [ ! -e ${dwimask}  ]; then
        echo "dwi2mask: ${dwimask}"
        dwi2mask ${dwi} ${dwimask}
        echo "    dwi2mask complete"
fi

###### Generate FODs using single-shell Constrained Spherical Deconvolution
if [ ! -e ${wmfods}  ]; then
    ###### Generate Response Function
    if [ ! -e ${rfvoxels}  ]; then
        echo "dwi2response"
        dwi2response ${rf_algo} ${dwi} ${rf} 
        echo "    dwi2response complete"
    fi
    echo "dwi2fod"
    dwi2fod csd -mask ${dwimask} ${dwi} ${rf} ${wmfods}
    echo "    dwi2fod complete"
fi


###### Generate tracks using WM_FODs
if [ ! -e ${tck}  ]; then
    echo "tckgen (${PREFIXFOD}): ${tck} using ${wmfods}"
    tckgen ${wmfods} ${tck} -act ${fivett} -backtrack -crop_at_gmwmi -seed_dynamic ${wmfods} -maxlength 250 -select ${tckgen_num} -cutoff 0.06
fi
### If tck file doesn't exist, you probably ran out of memory
if [ ! -e ${tck}  ]; then
    echo "  tckgen did not finish"
fi
if [ -e ${tck}  ]; then
    echo "    tckgen completed successfully"
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    minutes=$(( ${DIFF} / 60 ))
    hours=$(( ${minutes} / 60 ))
    minutes=$(( ${minutes} - $(( 60 * ${hours} )) ))
    echo "****** Job execution has taken ${hours} hours and ${minutes} minutes so far"
fi

###### Reduce the number of tracks to save memory when SIFTing
if [ ! -e ${tck_reduced}  ]; then
    echo "tckedit (${PREFIXFOD}): ${tck} => ${tck_reduced}"
    tckedit -number ${tckedit_num} ${tck} ${tck_reduced}
fi
### If tck_reduced file doesn't exist, you probably ran out of memory
if [ ! -e ${tck_reduced}  ]; then
    echo "    tckedit did not finish"
fi
if [ -e ${tck_reduced}  ]; then
    echo "    tckedit completed successfully"
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    minutes=$(( ${DIFF} / 60 ))
    hours=$(( ${minutes} / 60 ))
    minutes=$(( ${minutes} - $(( 60 * ${hours} )) ))
    echo "****** Job execution has taken ${hours} hours and ${minutes} minutes so far"
fi

###### SIFT tck_reduced file
if [ ! -e ${sft}  ]; then
    echo "tcksift (${PREFIXFOD}): ${tck_reduced} => ${sft}"
    tcksift -act ${fivett} -term_number ${tcksift_num} ${tck_reduced} ${wmfods_sft} ${sft}
fi
### If sft file doesn't exist, you almost certainly ran out of memory
if [ ! -e ${sft}  ]; then
    echo "  tcksift did not finish"
fi
if [ -e ${sft}  ]; then
    echo "    tcksift completed successfully"
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    minutes=$(( ${DIFF} / 60 ))
    hours=$(( ${minutes} / 60 ))
    minutes=$(( ${minutes} - $(( 60 * ${hours} )) ))
    echo "****** Job execution has taken ${hours} hours and ${minutes} minutes so far"
fi

###### HOORAY! ######
echo "****** Job completed (see above for errors) ******"
date
END=$(date +%s)
DIFF=$(( $END - $START ))
minutes=$(( ${DIFF} / 60 ))
hours=$(( ${minutes} / 60 ))
minutes=$(( ${minutes} - $(( 60 * ${hours} )) ))
echo "****** Job execution took ${hours} hours and ${minutes} minutes"
    """ % (h_vmem, tmem, nImages, jobName, workingDir, '%s', algorithm, workingDir, textFileListOfFilenames, '%s', '%s', '%s', '%s')
    # .format(h_vmem, tmem, nImages, workingDir, textFileListOfFilenames)
    qsubFile = textFileListOfFilenames + '.sh'
    outfile = open(qsubFile, 'w')
    print(qsubText,file=outfile)
    outfile.close()
    return qsubFile


def main():
    """
    Put it all together and hope for the best. ;-)
    
    Modify as desired: 
      algorithm = gif2 or gif3 or freesurfer
      workingDir
    """
    #Tiantian: change the dir to yours
    workingDir = '/cluster/project9/HCP_subjects/anatomically_constrained_tractography/ADNI/act'
    folderContaining5TTandWMFODs = sys.argv[1]
    JOBNAME = folderContaining5TTandWMFODs.replace('_','')
    algorithm = sys.argv[2] #gif3
    
    #algorithm = 'gif2'
    #JOBNAME = "HCPGIF2SC"
    #algorithm = 'freesurfer'
    #JOBNAME = "HCPFSSC"
    
    topFolder = os.path.normpath(os.path.join(os.path.normpath(workingDir),folderContaining5TTandWMFODs))
    
    runDate = str(datetime.date.today()).replace('-','')
    #*** Raw images 
    wildcardRaw = '*_5tt_{0}.mif'.format(algorithm)
    #wildcardRaw = '*_WMFODs_tournier.mif'.format(algorithm)
    depth = 1 # use glob in the top level directory only
    filez = traverseFolders(topFolder,wildcardRaw,depth)
    filez.sort()
    filez_prefix = traverseFolders(topFolder,wildcardRaw,depth)
    #filez_prefix = ['_'.join(os.path.basename(pathToImage).split('_')[0:3]) for pathToImage in filez]
    filez_prefix = [s.replace(wildcardRaw.replace('*',''),'') for s in filez_prefix]
    print(filez_prefix)
    #*** Processed images
    wildcardProcessed = '*tcksift*.tck'
    depth = 1 # anything other than 1 uses os.walk down through the directory tree
    filezDone = traverseFolders(topFolder,wildcardProcessed,depth)
    filezDone.sort()
    filezDone_prefix = ['_'.join(os.path.basename(pathToImage).split('_')[0:3]) for pathToImage in filezDone]
    filezDone_prefix5TT = [os.path.basename(s).replace(wildcardRaw,'') for s in filezDone_prefix]
    filezDone_prefixFOD = [os.path.basename(s)[5:15] for s in filezDone_prefix]
    
    print("= = = = = = = = = = = = filezDone_prefix: (move these files to Done subfolder) = = = = = = = = = = = =")
    print(filezDone_prefix5TT)
    print(filezDone_prefixFOD)
    
    #*** Move processed files
    dest = '{0}/Done'.format(topFolder)
    print('Moving processed files to Done/ subfolder')
    for k in range(0,len(filezDone_prefix)):
        filesToMove_tck = '{0}/{1}*.tck'.format(topFolder,filezDone_prefix[k])
        filesToMove_mif = '{0}/{1}*.mif'.format(topFolder,filezDone_prefix[k])
        #tx = os.system('mv {0} {1}/'.format(filesToMove_tck,os.path.normpath(dest)))
        #tx = os.system('mv {0} {1}/'.format(filesToMove_mif,os.path.normpath(dest)))
        #print 'Moved {0} and {1}\n'.format(filesToMove_tck,filesToMove_mif)
        print('=== Move {0} and {1}\n'.format(filesToMove_tck,filesToMove_mif))
    
    #*** Yet to be processed
    #filezYetToBeProcessed = [pathToImage for pathToImage in filez_prefix if pathToImage not in filezDone_prefix]
    ##filezYetToBeProcessed = [pathToImage + '_5TT_' + algorithm + '.mif' for pathToImage in filezYetToBeProcessed]
    #filezYetToBeProcessed = [pathToImage + '_WMFODs_tournier.mif' for pathToImage in filezYetToBeProcessed]
    filezYetToBeProcessed = filez
    
    # Write list of image files for processing to text file
    filenameList = '{0}/{1}_filelist_{2}.txt'.format(topFolder,JOBNAME,runDate)
    outfile = open(filenameList, 'w')
    print("\n".join(str(i) for i in filezYetToBeProcessed),file=outfile)
    outfile.close()
    
    # Generate qsub file
    qsubFile = qsubSC(topFolder,filenameList,JOBNAME,algorithm)
    # Print instructions for submitting the SC pre-processing to the cluster
    print("\n * * * Found {0} subjects:\n       - {1} processed (based on filename containing tcksift)\n       - {2} prepared for cluster submission\n * * * \n".format(len(filez),len(filezDone),len(filezYetToBeProcessed)))
    print('Submit your job using:\n  qsub {0}'.format(qsubFile))
    
    # Write list of processed image files to text file
    filenameList = filenameList.replace('.txt','_done.txt')
    outfile = open(filenameList, 'w')
    print("\n".join(str(i) for i in filezDone),file=outfile)
    outfile.close()


if __name__ == '__main__':
    main()


