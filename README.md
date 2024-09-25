[![DOI](https://zenodo.org/badge/260175954.svg)](https://zenodo.org/badge/latestdoi/260175954)

# mcera5 <img src="inst/figures/hex.png" align="right" height="200"/>

A package to download and process ERA5 data ready for use in microclimate modelling. Corresponding paper describing the package [here in _Methods in Ecology and Evolution_](https://doi.org/10.1111/2041-210X.13877).

## Install

You can install the package from this GitHub repository via the follow line:  
`remotes::install_github("dklinges9/mcera5")`

## News

_Sept 2024_: This package is going through a lot of development to keep up with changes as ECMWF migrates to the new [beta Climate Data Store](https://cds-beta.climate.copernicus.eu/). These changes include lower API limits, different variable names/netCDF file structures, and different credentials (see below for details). If you are facing unexplainable errors, please bring them to my attention (see Questions, Concerns, Issues below).  

_25 Sept 2024_: ERA5-reanalysis and ERA5-land files downloaded from the beta CDS have a different structure to their time dimension. To extract climate data from these files (e.g. using `extract_clim()`), please see the developer branch "new_cds_time_edits". You can install this branch directly via:

`remotes::install_github("dklinges9/mcera5", ref = "new_cds_time_edits")`

## Tutorial

The vignette can then be accessed via:   
`vignette("mcera5_vignette")`  

Alternatively you can navigate to the vignette manually [here](https://github.com/dklinges9/mcera5/blob/master/vignettes/mcera5_vignette.Rmd) on the GitHub repository.  

## Setup

ERA5 climate data can be downloaded from the ECMWF climate data store (CDS). Note that in July 2024 the CDS migrated to a new platform, and the old platform was deprecated in Sept 2024. The following describes how to access the data using R:

1) Register for an ECMWF account [here](https://accounts.ecmwf.int/auth/realms/ecmwf/login-actions/registration?client_id=cds&tab_id=IA1LKqgLVc0). Upon registering you will need to accept all of the Terms and Conditions listed at the bottom of the form.

2) Then, navigate to the CDS site [here](https://cds-beta.climate.copernicus.eu/) and login using the button in the top right. Once logged in, hover your mouse over your name in the top right, and click on the option "Your profile" that appears (this should bring you to [this page](https://cds-beta.climate.copernicus.eu/profile). Here you will 
find your User ID (UID) and Personal Access Token, both which are required for you to remotely download data from the CDS. Make a note of these.  

3) Each CDS dataset has its own unique Terms of Use. You will need to accept these Terms for ERA5-reanalysis at [this page](https://cds-beta.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download) (scroll down to "Terms of use" and accept). This same set of terms also applies for other Copernicus products, including ERA5-land.

## Questions, Concerns, Issues

Before emailing about a concern, please submit a reproducible example as an issue on this Github repository. To do so, navigate to [the main page of this repository](https://github.com/dklinges9/mcera5/) and click on the "Issues" tab, and then "New Issue." If you do not have a GitHub account or are otherwise unable to submit an issue, then contact David Klinges at dklinges9@gmail.com.

## Contact and Contributors

David Klinges: _Maintainer, primary contact_ (dklinges9@gmail.com)  
James Duffy: _Creator_  
Ilya Maclean: _Contributor_  
