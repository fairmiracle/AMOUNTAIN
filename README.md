Bioconductor page http://bioconductor.org/packages/AMOUNTAIN

# Installation
Install the stable version from Bioconductor
```
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("AMOUNTAIN")
```

Install the developer version from github, on Linux with GSL:
```
git clone https://github.com/fairmiracle/AMOUNTAIN.git
cd AMOUNTAIN
gcc -c src/AMOUNTAIN.c -fPIC -std=c99
gcc -shared -o src/AMOUNTAIN.so AMOUNTAIN.o -lgsl -lgslcblas
rm AMOUNTAIN.o
```

To use C version functions in R:
```
source('R/AMOUNTAIN.R')
dyn.load(paste("src/AMOUNTAIN", .Platform$dynlib.ext, sep = ""))
source('R/AMOUNTAINC.R')
```
Here is a table of C-version functions and pure R functions:

| C-version|      Pure R   |  Brief description                  |
|:----------|:-------------|:--------------------------------------------|
| `CGPFixSS` |  `moduleIdentificationGPFixSS` | Module identification on single network |
| `CGPFixSSTwolayer` | `moduleIdentificationGPFixSSTwolayer` | Module identification on two-layer network |
| `CGPFixSSMultiLayer` |  `moduleIdentificationGPFixSSMultilayer` | Module identification on multi-layer network |

# Refererence
Dong Li, Shan He, Zhisong Pan, Guyu Hu. Active modules for multilayer weighted gene co-expression networks: 
a convex optimization approach. biorxiv 2016. http://www.biorxiv.org/content/early/2016/06/03/056952
