# Analysis of sea Storm near Liverpool

## Structure of files

Directories: 
- Data.New:        Dada generated after analysis.
- Data.Processed:  Data of ocean climate in csv after some processing. 
- Data.Raw:        Data of ocean climate as probided.
- Figures:         Figures generated by the different analyses. 
- Report:          Source files for the file Report.pdf
- R.functions:     Several files with R code for auxiliary functions.
- R.old.code:      Old R scripts that are not used in the final analyses. 

Files: 
- Analysis_Clusters.R: analyses storm clustering before and after 
     applying harmonic analysis to remove periodicities. 
- Analysis_Storms.R: A sample script where the storm characteristics 
     from Callaghan et. al. are extracted from the Data.New files of storms. 
- Clean_Data.R: reads files from Data.Raw and rewrites them into Data.Processed 
     in a uniform csv format. 
- Modelling_Clusters.R: analysis of the different modelling options for the 
     interarrival time between storms.
- Generate_Storms.R: program that sample of future storms acording to 
      some different scenarios. 
- Report.pdf : Pdf file in article like format to describe the different modelling 
     options, performance and final choices made. 


- Compare-Harm-Poiss.R: program to compare Harm-Analisis.R and PoisProcAnalisis.R outputs. 
- PoisProcAnalisis.R: Analisis of interarrival time of storms via a non-homogenuous Poisson Process. 


TODO: put the scripts in a way that they really do what is said in this files. Important reminders: 

Analysis_Storms: figures. -> year date table of events before and after 
                        harmonic analysis. Index of variation vs threshold? 
Modelling_Clusters: Compare 4 methods 1. raw homogenous Poisson, 2. non-homogenuous Poisson, 3 GAM models, 4 Stationary Density kernel of the empirical distribution after 

Generate_Storms: generate storms under by the best models in the 
modelling part (3 different samples), but also analyse sensitivity
with respect to climate chage and clustering. 
1. Climate change alternative scenarios Scenarios: a: rpc26, b: rpc85. 
2. Alternative clustering scenarios: 
     a. storms occur equispaced in time at same rate of storms/year
     b. storms occur modelled spaced in time as a poisson distribution 
              at the same rate of stroms/year. 

Out of sample error in modelling of clusters via checking the storms in the dates where 
the tide is not available. 


