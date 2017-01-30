[![platform](http://www.bioconductor.org/shields/availability/devel/AMOUNTAIN.svg)](https://www.bioconductor.org/packages/devel/bioc/html/AMOUNTAIN.html#archives) [![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/AMOUNTAIN.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/AMOUNTAIN/)

Bioconductor page http://bioconductor.org/packages/AMOUNTAIN

# Installation
Install the stable version from Bioconductor
```
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("AMOUNTAIN")
```

Install the development version from Github
```
library(devtools)
install_github("fairmiracle/AMOUNTAIN")
```
## Compile on Linux
Make sure GSL in installed, type `gsl-config` in terminal.

To compile C code:
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

## Compile on Windows
It is not that straightforward to compile with GSL under Windows. Someone has created GSL Windows DLL and headers for both 32 and 64-bit in https://code.google.com/archive/p/oscats/downloads. Extract gsl-1.15-dev-win32.zip and gsl-1.15-dev-win64.zip into two directories:

 - C:\GSL\i386
 - C:\GSL\x64

and set the environment variable `LIB_GSL` as  `C:/GSL` instead of `C:\GSL`. Finally add `C:\GSL\x64\bin` to the `Path` in case missing the DLLs. Then the source can be compiled.

# Refererence
Dong Li, Shan He, Zhisong Pan, Guyu Hu. Active modules for multilayer weighted gene co-expression networks: a convex optimization approach. biorxiv 2016. http://www.biorxiv.org/content/early/2016/06/03/056952
