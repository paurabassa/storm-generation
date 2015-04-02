# Analysis of sea Storm near Liverpool

## Structure of files

Directories: 
- Raw.Data:    Data of ocean climate as probided
- Clean.Data:  Data of ocean climate in csv after some processing. 
- New.Data:    Dada generated after analysis
- R.functions: Several files with R code of auxiliary functions 

Files: 
- clean_data.R: reads files from Raw.Data and rewrites it into Clean.Data in
                uniform csv format. 
- analysis_clusters.R: does some plots to analyise the clustering of storm events. 
- Compare-Harm-Poiss.R: program to compare Harm-Analisis.R and PoisProcAnalisis.R outputs. 
- Harm-Analisis.R: Analisis of interarrival time of storms via removing periodic components in hs
- PoisProcAnalisis.R: Analisis of interarrival time of storms via a non-homogenuous Poisson Process. 

