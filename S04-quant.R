## Quantitative proteomics

##      | Label-free | Labelled    |
## -----+------------+-------------|
## MS1  | XIC        | SILAC, 15N  |
## MS2  | Counting   | TMT (iTRAQ) |

## https://lazear.github.io/sage/
## https://uclouvain-cbio.github.io/sager/

library(tidyverse)
library(msdata)
library(QFeatures)

## Quantitative assay -> SummarizedExperiment
##
## - quantitative data: matrix [features x samples]
## - samples/columns annotation (meta-data)
## - features/rows annotation (meta-data)

## - quantitative data: assay()
## - samples/columns annotation: colData()
## - features/rows annotation: rowData()
