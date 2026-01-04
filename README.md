# code_Why-sensory-neurons-are-tuned-to-multiple-stimulus-features

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18143325.svg)](https://zenodo.org/records/18143325)

## Environment

* **Matlab 2024a**

## Data Availability

The datasets required for this project are hosted on Zenodo due to file size limits.
**Download Link:** [https://zenodo.org/records/18143325](https://zenodo.org/records/18143325)

> **Note:** Please download the data and ensure the `.mat` files are placed in their corresponding folders as described below to ensure the scripts run correctly.

## Data Preparation & Usage

Below is the guide to reproducing the specific figures using the scripts provided in each folder.

### 1. Figure 1B(right) & Figure 1C(right)
* **Folder:** `francsyn`
* **Script to Run:** `info_calculation_forL.m`
* **Data Dependencies:**
    The script reads the following data (downloaded from Zenodo) from the `francsyn` directory:
    * `neuron_data.mat`
    * `toyvscell_fullresponse.mat`

### 2. Figure 3C & Figure 3D
* **Folder:** `behaviral_synergy`
* **Script to Run:** `behav_inf.m`
* **Data Dependencies:**
    The script reads the following data from the `behaviral_synergy` directory:
    * `ga_allinfodata.mat`
    * `re allinfodata.mat`
    * `xt allinfodata.mat`

### 3. Figure 4C-E
* **Folder:** `opulation Stimulus synergy`
* **Script to Run:** `MT_popss.m`

### 4. Figure 5A-D
* **Folder:** `2D_Behavior_Synergy`
* **Scripts to Run:**
    * `popss5ab.m`
    * `popss5cd.m`

### 5. Figure 6A-C
* **Folder:** `MT_popsim_multiD`
* **Script to Run:** `MT_popsim_multiD_encoding.m`
