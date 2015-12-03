[![Build Status](https://travis-ci.org/pblischak/polyfreqs.svg?branch=master)](https://travis-ci.org/pblischak/polyfreqs)

# **polyfreqs**

## An R package for Bayesian population genomics in autopolyploids

**polyfreqs** is an R package for the estimation of biallelic SNP frequencies, genotypes and heterozygosity in autopolyploid taxa using high throughput sequencing data.


### Dependencies

**polyfreqs** uses C++ code to implement its Gibbs sampling algorithm which will usually require the installation of additional software (depending on the operating system [OS] being used).
Windows users will need to install <a href="http://cran.r-project.org/bin/windows/Rtools/" target="_blank">Rtools</a>.
MacOSX users will need to install the Xcode Command Line Tools.
Linux users will need an up-to-date version of the GNU C Compiler (gcc) and the r-base-dev package.  
Looking at the requirements for Rcpp is a good place to start too.

**polyfreqs** relies on two other R packages: <a href="http://cran.r-project.org/package=Rcpp" target="_blank"><strong>Rcpp</strong></a> and <a href="http://cran.r-project.org/package=RcppArmadillo" target="_blank"><strong>RcppArmadillo</strong></a>.
These are both available on CRAN and can be installed in the usual way using `install.packages()`:

```r
install.packages("Rcpp")
install.packages("RcppArmadillo")
```

Note that **Rcpp** and **RcppArmadillo** also require the compilation of C++ code so make sure that the necessary compilers are installed appropriately for your OS.

### Installation

**polyfreqs** v1.0.0 is now on CRAN: <a href="http://cran.r-project.org/package=polyfreqs" target="_blank">link</a>.

You can now install it like you would any other R package:

```r
install.packages("polyfreqs")
```

Installing the latest developmental release of **polyfreqs** can be done using the <a href="http://cran.r-project.org/package=devtools" target="_blank"><strong>devtools</strong></a> package and the `install_github()` command.
Install **devtools** using `install.packages("devtools")`. **polyfreqs** can then be installed as follows:

```r
devtools::install_github("pblischak/polyfreqs")
```

### Documentation

Example code and tutorials for running **polyfreqs** can be found in the <a href="https://cran.r-project.org/web/packages/polyfreqs/vignettes/polyfreqs_Intro.html" target="_blank">vignette</a>. 
For more details on the model underlying **polyfreqs** please see the associated paper in *Molecular Ecology Resources*: <a href="http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12493/abstract" target="_blank">Blischak <em>et al</em>.</a> The Supplemental Material also has a walk through for analysis a data set collected for autotetraploid potato (*Solanum tuberosum*).

--------

**Release notes**

 - **v1.0.0** -- First release. Now available on CRAN.
