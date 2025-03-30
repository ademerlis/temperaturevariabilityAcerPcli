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
* **.Rmd:** has all the code needed for Figure 4 and survivorship analysis. It sources te data from the **Data** subfolder and saves statistical summaries and plots in the **Outputs** subfolder.

#### 4.YII:
* **1.IPAM_ImportFunction.R:** R function to imports the csv raw IPAM data. Yo do not have to do anything with this file, but the function is called by this file below.
* **2.MergeYIIandID.R:** Applies **1.IPAM_ImportFunction.R** to the raw data in the **IPAM_Raw**, adds the sample metadata in the **ID_AOI.csv** file and creates a long format file with the YII values that is exported to the file **Data/YII_tall.csv**
* **3.YII_Acer_Nut.Rmd:** has all the code needed for Figure 2 and related statistical analisis. It sources the data from the file **Data/YII_tall.csv** and saves statistical summaries and plots in the **Outputs** subfolder.

#### 5.Symbiodiniaceae:
* **SH_cell_code.Rmd:** has all the code needed for Figure 3. It sources the qPCR data from the **Data** subfolder and the sample metadata from the **Metadata.csv** file. Statistical summaries and plots are saved in the **Outputs** subfolder.
* **Data:** has all the raw data exported from the qPCR machine.

#### 6.Metadata:
* This folder has copies of the metadata files, but are actually not used by any code. 

</br>

#### Notes
* This README.md was formatted following [Dr. Ana Palacio's GitHub repository](https://github.com/anampc/Acer_NH4_disease/tree/master). 
