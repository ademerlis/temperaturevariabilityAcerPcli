#[![DOI](https://zenodo.org/.svg)](https://zenodo.org/doi)

This repository contains data and analysis scripts for the manuscript:

## Species-level differences in molecular responses to a thermally variable stress-hardening treatment for Caribbean corals
#### **Authors:** Allyson DeMerlis, Michael S. Studivan, Kevin Wong, Nash Soderberg, David Ehrens, Lys M. Isma, Katrina Rosing, Katrina Sophia Cocson, Rowan Thomas, Danielle Dvorkin, Patrick M. Kiel, Joseph D. Unsworth, Martine D’Alessandro, Ana M. Palacio-Castro, Diego Lirman, Andrew C. Baker, Erinn M. Muller, Nikki Traylor-Knowles, Ian C. Enochs
#### **Journal:** _Submitted to Molecular Ecology_ [doi:XXX](http://dx.doi.org/XXX)  

-----

### Description:
These repository contains all data and code used to study the impact of thermal variability as an ex-situ stress-hardening treatment on the physiology and gene expression of _Acropora cervicornis_ and _Pseudodiploria clivosa_.

### Contents:

#### Tank Parameters:
* **tank_temp_figs.Rmd:** has all the code needed for Figure 1. It sources the data from the **Data from LabVIEW** subfolder and saves statistical summaries and plots in the **plots** subfolder.

#### Physiology:
* The complete data file that contains all metadata and buoyant weight mesaurements is the file **metadata.csv**.
* Physiology data analysis is broken up into three subfolders: **Calcification**, **R_intensity**, and **Photosynthetic Efficiency**.
* In the subfolder **Calcification**, **calcification.Rmd** has all the code needed for Supplementary Figure 1 and Table S3 statistics.
* In the subfolder **R_intensity**, **colorscores.Rmd** has all the code needed for Supplementary Figure 2 and Table S4 statistics.
* The subfolder **Photosynthetic Efficiency** has several R files that are important for analysis.
  - **1_importtidyPAMdata.Rmd** includes the function to import raw IPAM data into R, and creates a format file that matches coral fragment metadata to the area of interest (AOI) and YII (photosynthetic efficiency) values to the IPAM image metadata files. The raw files for this can be found in the **ipam_data** subfolder.
  - 


#### 5.Symbiodiniaceae:
* **SH_cell_code.Rmd:** has all the code needed for Figure 3. It sources the qPCR data from the **Data** subfolder and the sample metadata from the **Metadata.csv** file. Statistical summaries and plots are saved in the **Outputs** subfolder.
* **Data:** has all the raw data exported from the qPCR machine.

#### 6.Metadata:
* This folder has copies of the metadata files, but are actually not used by any code. 

</br>

#### Notes
* This README.md was formatted following [Dr. Ana Palacio's GitHub repository](https://github.com/anampc/Acer_NH4_disease/tree/master).

