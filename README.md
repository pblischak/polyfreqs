# **polyfreqs**

## Estimating allele frequencies in populations of autopolyploids

**polyfreqs** is an R package for the estimation of biallelic SNP frequencies in polyploid taxa. 


### Dependencies

**polyfreqs** uses C++ code to implement its Gibbs sampling algorithm which will usually require the installation of additional software (depending on the operating system [OS] being used). 
Windows users will need to install [Rtools](http://cran.r-project.org/bin/windows/Rtools/). 
MacOSX users will need to install the Xcode Command Line Tools. 
Linux users will need an up-to-date version of the GNU C Compiler (gcc) and the r-base-dev package.  
Looking at the requirements for Rcpp is a good place to start too.

**polyfreqs** relies on two other R packages: [**Rcpp**](http://cran.r-project.org/web/packages/Rcpp/index.html) and [**RcppArmadillo**](http://cran.r-project.org/web/packages/RcppArmadillo/index.html). 
These are both available on CRAN and can be installed in the usual way using `install.packages()`:

```r
install.packages("Rcpp")
install.packages("RcppArmadillo")
```

Note that **Rcpp** and **RcppArmadillo** also require the compilation of C++ code so make sure that the necessary compilers are installed appropriately for your OS.

### Installation

Installing **polyfreqs** can be done using the [**devtools**](http://cran.r-project.org/web/packages/devtools/index.html) package and the `install_github()` command. 
Install **devtools** using `install.packages("devtools")`. **polyfreqs** can then be installed as follows:

```r
devtools::install_github("pblischak/polyfreqs")
```

### Documentation

Example code and tutorials for running **polyfreqs** can be found in the <a href="https://github.com/pblischak/polyfreqs/wiki" target="_blank">wiki</a>. 
For more details on the model underlying **polyfreqs** please see the associated manuscript on bioRxiv: <a href="https://wolfelab.wordpress.com" target="_blank">Blischak <em>et al<\em>. 2015</a>.
