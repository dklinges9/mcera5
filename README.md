[![DOI](https://zenodo.org/badge/260175954.svg)](https://zenodo.org/badge/latestdoi/260175954)

# mcera5 <img src="inst/figures/hex.png" align="right" height="200"/>

A package to download and process ERA5 data ready for use in microclimate modelling. Corresponding paper describing the package [here in _Methods in Ecology and Evolution_](https://doi.org/10.1111/2041-210X.13877).

## Install

At present, `mcera5` has a dependency on the R package `microclima`, which in turn depends on `rgdal`. Given that `rgdal` has been retired, it is no longer available to be installed from CRAN. Thus, if you do not have `rgdal` available for your current version of R, _prior_ to attempting to install `mcera5`, please install `rgdal` from an archived source by running this line of code:  
`remotes::install_url("https://cran.r-project.org/src/contrib/Archive/rgdal/rgdal_1.6-7.tar.gz", type="source")`

Then when installing `mcera5`, make sure to specify `build_vignettes = TRUE` in order to access the package's vignette tutorial:  
`remotes::install_github("dklinges9/mcera5", build_vignettes = TRUE)`

## Tutorial

The vignette can then be accessed via:   
`vignette("mcera5_vignette")`  

Alternatively you can navigate to the vignette manually [here](https://github.com/dklinges9/mcera5/blob/master/vignettes/mcera5_vignette.Rmd) on the GitHub repository.  

## Contact and Contributors

David Klinges: _Maintainer, primary contact_ (dklinges9@gmail.com)  
James Duffy: _Creator_  
Ilya Maclean: _Contributor_  

## Issues

Please submit reproducible examples as an issue on Github. 
