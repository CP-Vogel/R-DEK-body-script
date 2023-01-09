# R-DEK-body-script
R script to aggregate time-resolved data from a microscopic time-lapse experiment

## Prerequisites

This script was developed with R software version 3.4.3

Following R libraries are required to run the script:

```
library(plyr)          
library(tidyverse) 
library(Hmisc)  
library(broom)     
library(data.table)
library(mmand)    
library(dtplyr) 
```

This script only accepts CSV files created by the CellProfiler image analysis pipeline (version 2.1.0)
