# Anatomically Constrained Tractography for Neurodegenerative Disease Progression Modelling

[Neil Oxtoby](https://github.com/noxtoby), UCL Centre for Medical Image Computing
February 2022

## History
This code was used to generate healthy elderly connectomes for network spreading models of pathology. 

It was develop using [PPMI](https://ppmi-info.org) data and the [MRTrix3 ACT pipeline](https://mrtrix.readthedocs.io/en/dev/quantitative_structural_connectivity/act.html) (see also [this ISMRM tutorial](https://mrtrix.readthedocs.io/en/dev/quantitative_structural_connectivity/ismrm_hcp_tutorial.html) on [HCP](http://www.humanconnectomeproject.org/) data).


## Contents

- `PPMI` folder: for data, including outputs
- `scripts` folder: the code
  - `act_ppmi_0_prep.py`: data wrangling
  - `act_ppmi_1_preprocess.py`: preprocessing steps to run locally
  - `ppmi_tckgen_cluster.py`: tractography to run on the UCL CMIC/CS cluster
  - `act_utilities.py`: utilities

PPMI imaging data was downloaded as raw ("Original") DICOM files. Will require tweaking for other datasets, such as ADNI.
