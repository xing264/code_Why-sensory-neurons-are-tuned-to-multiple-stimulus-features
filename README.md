# code_Why-sensory-neurons-are-tuned-to-multiple-stimulus-features

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18143325.svg)](https://zenodo.org/records/18143325)

## Environment

* **Matlab 2024a**

## Data Availability

The data for this project is hosted in two locations depending on file size:

1.  **Included in Repository:** Small datasets (e.g., for Figure 1) are directly included in the corresponding folders in this GitHub repository.
2.  **Hosted on Zenodo:** Larger datasets (e.g., for Figure 3/Behavioral Synergy) are hosted on Zenodo.
    * **Download Link:** [https://zenodo.org/records/18143325](https://zenodo.org/records/18143325)

> **Note:** For data downloaded from Zenodo, please place the `.mat` files into their corresponding local folders as indicated below.

## Data Preparation & Usage

Below is the guide to reproducing the specific figures using the scripts provided in each folder.

### 1. Figure 1B(right) & Figure 1C(right)
* **Folder:** `francsyn`
* **Script to Run:** `info_calculation_forL.m`
* **Data Source:** **Included in this repository.**
    * The script reads the local files: `neuron_data.mat` and `toyvscell_fullresponse.mat`.

### 2. Figure 3C & Figure 3D
* **Folder:** `behaviral_synergy`
* **Script to Run:** `behav_inf.m`
* **Data Source:** **Download from Zenodo.**
    * Please download and place the following files into the `behaviral_synergy` directory:
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
