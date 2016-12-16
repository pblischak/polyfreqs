[![Build Status](https://travis-ci.org/pblischak/polyfreqs.svg?branch=master)](https://travis-ci.org/pblischak/polyfreqs)
[![](http://www.r-pkg.org/badges/version/polyfreqs)](https://CRAN.R-project.org/package=polyfreqs)
[![CRAN downloads](http://cranlogs.r-pkg.org/badges/polyfreqs)](https://CRAN.R-project.org/package=polyfreqs/vignettes/polyfreqs_Intro.html)

# **polyfreqs**

## An R package for Bayesian population genomics in autopolyploids

**polyfreqs** is an R package for the estimation of biallelic SNP frequencies, genotypes and heterozygosity in autopolyploid taxa using high throughput sequencing data. It should work for diploids as well, but does not accomodate data sets of mixed ploidy.

### **polyfreqs** does accept missing data: code them as `0` in the total read count matrix.

> **NEW**: **polyfreqs** now has a Google Groups page. Please feel free to join the group and post any questions that you may have about the software. [[Google Groups link](https://groups.google.com/forum/#!forum/polyfreqs-users)]

### Dependencies

**polyfreqs** uses C++ code to implement its Gibbs sampling algorithm which will usually require the installation of additional software (depending on the operating system [OS] being used).
Windows users will need to install <a href="https://CRAN.R-project.org/bin/windows/Rtools/" target="_blank">Rtools</a>.
MacOSX users will need to install the Xcode Command Line Tools.
Linux users will need an up-to-date version of the GNU Compiler Collection (gcc) and the r-base-dev package. **polyfreqs** relies on the R package <a href="https://CRAN.R-project.org/package=Rcpp" target="_blank"><strong>Rcpp</strong></a> which is a good place to start too for figuring what you will need. Note that **Rcpp** also requires the compilation of C++ code so make sure that the necessary compilers are installed appropriately for your OS. You can install **Rcpp** directly from CRAN in the usual way using the `install.packages()` command:

```r
install.packages("Rcpp")
```

### Installation

**polyfreqs** v1.0.0 is now on CRAN: <a href="https://CRAN.R-project.org/package=polyfreqs" target="_blank">link</a>.

You can now install it like you would any other R package:

```r
install.packages("polyfreqs")
```

Installing the latest developmental release of **polyfreqs** can be done using the <a href="https://CRAN.R-project.org/package=devtools" target="_blank"><strong>devtools</strong></a> package and the `install_github()` command.
Install **devtools** using `install.packages("devtools")`. **polyfreqs** can then be installed as follows:

```r
devtools::install_github("pblischak/polyfreqs")
```

### Documentation

Example code and tutorials for running **polyfreqs** can be found in the <a href="https://CRAN.R-project.org/package=polyfreqs/vignettes/polyfreqs_Intro.html" target="_blank">vignette</a>.
For more details on the model underlying **polyfreqs** please see the associated paper in *Molecular Ecology Resources*: <a href="http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12493/abstract" target="_blank">Blischak <em>et al</em>.</a> The Supplemental Material also has a walk through for analyzing a data set collected for autotetraploid potato (*Solanum tuberosum*).

--------

**Release notes**

 - **v1.0.2** -- Small patch that updated code for sampling genotypes during the MCMC that was giving underflow errors when total read counts are high (~1000x coverage).

 - **v1.0.1** -- Removed dependency on the **RcppArmadillo** `sample()` function by coding our own version (`nonunif_int()` in the **sample_g.cpp** source file). The Gibbs sampler should run a bit faster now.

 - **v1.0.0** -- First release. Now available on CRAN.
