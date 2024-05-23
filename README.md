# Novel machine learning model for predicting cancer drugs'susceptibilities and discovering novel treatments
---

## About this Project

When the relationship between predictors and the outcome is complex and nonlinear, multi-task prediction models based on linear assumptions tend to perform poorly. Kernel functions can project the data into a higher-dimensional space, thereby simplifying these complex nonlinear relationships. To address this issue, we developed an improved method based on MTPS, called Kernelized Residual Stacking (KRS).
---
## Directory Layout

We assume the user is running the script on R version **4.2.0** and has set up a directory with the same structure.
~~~
    [your_deirctory]  
~~~
R codes are in the subdirectory directory at **code**. 
~~~
    [your_deirctory]/code  
~~~
All numerical results saved during the code execution, which are used to produce graphs, are located in the **result**.
~~~
    [your_deirctory]/result 
~~~

all the graphs in the paper are in the subdirectory directory at **figure**. 
~~~
    [your_deirctory]/figure  
~~~
the supplementary figure in the paper is in the subdirectory directory at **supplementary figure**. 
~~~
    [your_deirctory]/supplementary figure  
~~~
all the **raw datasets** from CCLE and GDSC are in the subdirectory directory at **data**. 
~~~
    [your_deirctory]/data  
~~~

<details><summary>code</summary>

    ├── code  
    │   ├── 1-data_cleaning.R		    # clean the raw data 
    │   ├── 2-models.R		    # 5 different models
    │   ├── 3-enrichment analysis.R		    # enrichment analysis
    │   └── 4-plot.R		    # The graphs of the main file and additional files		
    
</details>
<details><summary>result</summary>

    ├── result
    │ 	 ├── csv files
    │ 	 │ 	 ├── ccle.csv		    # cmap analysis result from L1000 analysis platform 
    │ 	 │ 	 ├── gdsc.csv		    # cmap analysis result from L1000 analysis platform
    │ 	 │ 	 ├── xx1gene(up-down).csv		    # geneselected gene using to do Cmap analysis
    │ 	 │ 	 └── xx2gene(up-down).csv		    # geneselected gene using to do Cmap analysis 
    │ 	 └── rda files
    │ 	     ├── bias.abs.GDSCtrain.rda		    # mse and bias result of 5 models
    │ 	     ├── drug7_gene1067.rda		    # preprocessed data for modelling
    │ 	     ├── enrichmentplot.rda		    # result of the enrichment analysis from string website
    │ 	     └── selectedgene.rda		    # result of gene selection used to do enrichment analysis on the string website

</details>
<details><summary>figure</summary>

    ├──  figure  
    │ 	 ├── cmap.png
    │ 	 ├── enrichment analysis.png
    │ 	 ├── hubgene.png
    │ 	 └── train_GDSC_mse_figure.png	
    ├── supplementary figure
    │ 	 └── boxplot_train_GDSC_bias_krs and rs.png

</details>
<details><summary>data</summary>

    ├──  data
    │    ├── CCLE_DrugResponse.csv		    # drug response data for CCLE
    │    ├── CCLE_ExpData.csv		    # gene expression data for CCLE
    │    ├── GDSC_DrugResponse.csv		    # drug response data for GDSC
    │    └── GDSC_ExpData.csv		    # gene expression data for GDSC

</details>
<details><summary>updated MTPS package</summary>

    └── updated MTPS package		    # we provided the updated MTPS package
        ├── MTPS
        │    ├── DESCRIPTION
        │    ├── NAMESPACE
        │    ├── R
        │    │    ├── AUC.r
        │    │    ├── MTPS.R
        │    │    ├── checkMatch.R
        │    │    ├── createFolder.R
        │    │    ├── cv.MTPS.R
        │    │    ├── cvGlmnet2.R
        │    │    ├── cvMultiFit.R
        │    │    ├── list.learners.R
        │    │    ├── method.R
        │    │    ├── multiFit.R
        │    │    ├── residuals.R
        │    │    └── rsMultiFit.R
        │    ├── README.md
        │    ├── build
        │    │    └── vignette.rds
        │    ├── data
        │    │    ├── HIV.rda
        │    │    └── Internal.rda
        │    ├── inst
        │    │    ├── CITATION
        │    │    └── doc
        │    │        ├── Guide.R
        │    │        ├── Guide.Rmd
        │    │        └── Guide.html
        │    ├── man
        │    │    ├── AUC.Rd
        │    │    ├── HIV_data.Rd
        │    │    ├── Internal_data.Rd
        │    │    ├── MTPS-internal.Rd
        │    │    ├── cv.MTPS.Rd
        │    │    ├── list.learners.Rd
        │    │    ├── modify.parameter.Rd
        │    │    ├── multiFit.Rd
        │    │    ├── predict.MTPS.Rd
        │    │    ├── predict.multiFit.Rd
        │    │    └── revisedStacking.Rd
        │    └── vignettes
        │        └── Guide.Rmd
        └── README.md

</details>

---
## Notice

---
## Before you start

1. decide the path of [your_directory] to replicate our results;
2. create the subdirectories **code**, **data**, **figure**, **result**, **supplementary result** at [your_directory]；
3. allocate all relevant files into each subdirectory. Folders should look like the figure below:
![image](https://github.com/ubcxzhang/KRS/blob/main/illustration.png?raw=true)
4. Please run the files in the code folder in the following order: "1-data cleaning.R", "2-models.R", "3-enrichment analysis.R", "4-plot.R".
---

<details><summary>1-data cleaning.R</summary>
- Read the raw data from `[your_directory]/data`.
    
- Remove missing values.
  
- Filter for overlapping drugs and genes between the two datasets.
    
- The cleaned data was saved in `[your_directory]/result/rda files/drug7_gene1067.rda`.
  
</details>

<details><summary>2-models.R</summary>
- Using the cleaned data, `[your_directory]/result/rda files/drug7_gene1067.rda`, five multi-task prediction models were established: RS, KRS, L21, KMTrace, and KBMTL. 
    
- The mean squared error (MSE) of the predictions and the absolute difference in bias were stored in `[your_directory]/result/rda files/bias.abs.GDSCtrain.rda`.
  
</details>

<details><summary>3-enrichment analysis.R</summary>
- Perform gene selection on the cleaned data, `[your_directory]/result/rda files/drug7_gene1067.rda`. 
    
- The selected genes were stored in `[your_directory]/result/csv files/xx1gene(up-down).csv` and `[your_directory]/result/csv files/xx2gene(up-down).csv`.
  
</details>

---
## Notice

Then, the GO analysis, KEGG analysis and PPI network were performed using **[your_directory]/result/csv files/xx1gene(up-down).csv** and **[your_directory]/result/csv files/xx2gene(up-down).csv** on String <https://string-db.org/>. Results were saved in **[your_directory]/result/rda files/enrichment plot.rda**. 

Cmap analysis was performed on L1000 analysis platform using the same csv files. Results were saved in **[your_directory]/result/csv files/ccle.csv** and **[your_directory]/result/csv files/gdsc.csv**.

---

<details><summary>4-plot.R</summary>
- Main figure of the MSE was created using `[your_directory]/result/rda files/bias.abs.GDSCtrain.rda`. 
    
- Main figure of the GO and KEGG were created using `[your_directory]/result/rda files/enrichmentplot.rda`.
  
- Main figure of PPI was created using the Cytoscape.
  
- Main figure of Cmap was created using `[your_directory]/result/csv files/ccle.csv` and `[your_directory]/result/csv files/gdsc.csv`.
  
- A supplementary figure for the difference between absolute bias was created using `[your_directory]/result/rda files/bias.abs.GDSCtrain.rda`.
  
       
</details>





