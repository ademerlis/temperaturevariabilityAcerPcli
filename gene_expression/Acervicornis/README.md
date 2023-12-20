# DESeq2 analysis for *Acropora cervicornis*

Last updated: 2023-12-20

## Results

Code to create these graphs is from [this R file](). 

`dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Treatment)`

Treatment has three levels:
1. Initial (day 0)
2. Untreated (control day 29)
3. Treated (variable treated day 29)


### Identifying outliers

Using the R package "biobase", I ran arrayQualityMetrics, which generates a report that visualizes which samples can be considered outliers.

![Screen Shot 2023-12-20 at 3 38 43 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/cbda63be-6894-4d2b-a2b4-abeebdcbbb90)


### PCoA plot

![Screen Shot 2023-12-20 at 3 58 10 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/23cb4d18-e89e-47c8-9c5d-54419c0de6fc)

### PERMANOVA

![Screen Shot 2023-12-20 at 4 01 12 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/dfcba9c5-c521-498b-bea3-22bacb3a8976)

### Significant DGEs

![Screen Shot 2023-12-20 at 4 12 29 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/6319c22e-ac59-41ea-a42e-da18eb06e4ef)

![Screen Shot 2023-12-20 at 4 12 59 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/aa7754a7-eb10-4518-a6d7-ae532304fc1c)

![Screen Shot 2023-12-20 at 4 14 19 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/08c140ab-1d8a-43ae-b658-ae7080981d4b)

![Screen Shot 2023-12-20 at 4 15 16 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/8282e08f-89d5-4ff3-a128-3da2f2f35095)

