# 2026 Tanner Crab CIE Review Repository
**NOTE: under construction**

This is the repository for materials for the 2026 CIE Review of the Bering Sea and Aleutian Islands (BSAI) Tanner crab assessment, which will be held June 9-11, 2026 at the Alaska Fisheries Science Center (Seattle, WA).
The BSAI Tanner crab (*Chionoecetes bairdi*) stock is one of [10 crab stocks](https://www.npfmc.org/fisheries/bsai-crab/) managed by the North Pacific Fishery Management Council ([NPFMC](https://www.npfmc.org/about-the-council/)). 
The review will include a number of topics including: 
  * assessment model transition from the current bespoke modeling framework (TCSAM02) to the framework used by the majority of 
BSAI crab stock assessments, the Generalized Model for Assessing Crustacean Stocks (GMACS)
  * issues related to management Tier level, management proxies for $F_{MSY}$ and $B_{MSY}$, and OFL and ABC calculations
  * issues related to whether, and if so how, hybrid *Chionoecetes* should be included in the assessment
  * issues related spatial distribution, stock structure, and management

## Review agenda

A pdf version of the current agenda is available at the following link:

  * [Agenda](documents/Tanner_Crab_CIE_Agenda_Table.pdf)

## Assessment-related documents

The 2025 Tanner crab assessment is available as a set of pdfs (main document, appendices, presentation) at: 

  * [2025 Tanner crab stock assessment](https://meetings.npfmc.org/Meeting/Details/3097)
  
The assessment describes the process by which federal management quantities and harvest levels were determined, as well as the then-current (Sept., 2025) state of 
transition from TCAM02 to GMACS.
  
The May 2026 Tanner crab report to the NPFMC's Crab Plan Team (CPT) consists of files at the following links:

  * [Report](documents/TannerCrabReport_2026-05.pdf)
  * [Addendum](documents/TannerCrabReport_2026-05_Addendum.pdf)
  * [Presentation](documents/TannerCrabReport_2026-05_Presentation.pdf)

The report and addendum discuss recent changes to GMACS that provide functionality equivalent to the bespoke Tanner crab modeling framework (TCSAM02) and the current state of 
matching models between the two frameworks. It also provides a preliminary discussion on the issue of including *Chionoecetes spp.* hybrids in the Tanner crab assessment.

### NMFS Eastern Bering Sea Trawl Survey

The NMFS eastern Bering Sea Trawl Survey (EBS trawl survey) is the primary non-fishery data source for the Tanner crab assessment. 
A comprehensive discussion of results from the 2025 survey can be found [here](documents/AFSC.EBS.CrabSurvey2025.NOAA-TM-AFSC_513.pdf).

## Assessment model frameworks

### TCSAM02 (the bespoke framework)

TCSAM02 is the bespoke size structured modeling framework for the Tanner crab assessment. Written in [ADMB](https://www.admb-project.org) tpl and C++, 
it was adopted in for the assessment 2017 and has been updated on at least an annual basis to 
address different issues and provide different options for the assessment. The code is hosted on [GitHub](https://github.com/wStockhausen/tcsam02/tree/202603); 
the current branch is '202603'. [wtsADMB](https://github.com/wStockhausen/wtsADMB) provides a library of ADMB-compatible C++ functions used in TCSAM02. 
A constellation of R packages (principally [rTCSAM02](https://github.com/wStockhausen/rTCSAM02) and [rCompTCMs](https://github.com/wStockhausen/rCompTCMs)) have 
evolved to extract and compare results between alternative models in various formats.

A detailed description of the TCSAM02 modeling framework is available [here](documents/TCSAM02Description2026.pdf).

Input files to run the TCSAM02 model are available in a zip file [here](model_files/TCSAM02/TCSAM02_InputFiles.zip). 
An OSX-compatible executable and shell script are available [here](model_files/TCSAM02/tcsam02) and [here](model_files/TCSAM02/run_tcsam02_osx.sh). 
A Windows-compatible executable and batch file are available [here](model_files/TCSAM02/tcsam02.exe) and [here](model_files/TCSAM02/run_tcsam02_bat.sh). 

To run the model (if you really want to):

  * extract files in the zip file to a folder (the "top-level folder")
  * copy the tcsam02 executable (tcsam02 or tcsam02.exe) to the top-level folder
    - you may have to set permissions to allow the executable to run
  * create a sub-folder "run" and copy the shell script or batch file into it
    - you may have to set permissions to allow the shell script/batch file to run
  * run the shell script/batch file from a terminal/command window opened in the sub-folder
  
The model should run twice: the first time to estimate model parameters at an approximate MLE then invert the hessian and estimate the covariance matrix, the second time to 
apply the "hess_step" procedure to determine the "true" MLE (all parameter gradients zero) then invert the hessian and estimate the covariance matrix and "std" file. 
Applying the "hess_step" procedure is probably overkill but provides peace-of-mind that the true (but still possibly local) maximum in the likelihood surface has been found.

If you have installed the necessary R packages, you can run the scripts in the "r_step1..." and "r_step2..." files to extract the results to an 
RData file and plot various quantities. The RData file (the result of running "r_step1...") can also be found [here](model_files/TCSAM02/Results.RData).

### GMACS (the "new" framework)

The Generalized Model for Assessing Crustacean Stocks, GMACS, is the modeling framework the assessment is in the process of transitioning to. It is also written in ADMB/C++. 
The principal reason for the transition is that GMACS is the modeling framework for the majority of BSAI crab stocks. It provides a single framework for both *Chionoecetes* and king crab life histories whereas TCSAM02 is specifically focused 
on Tanner crab. The [GMACS Project](https://github.com/GMACS-project) 
provides an organizing framework for repositories related to GMACS. The [GMACS_tpl-cpp_code](https://github.com/GMACS-project/GMACS_tpl-cpp_code) repository hosts 
the latest [version](https://github.com/GMACS-project/GMACS_tpl-cpp_code/tree/devel_202605) of GMACS, which includes several recent additions to provide features found 
in the bespoke model. While several R packages have been developed to extract and present results from GMACS models, 
the R package [wtsGMACS](https://github.com/wStockhausen/wtsGMACS) facilitates direct comparisons with TCSAM02 models.

Input files to run the GMACS model are available in a zip file [here](model_files/GMACS/GMACS_InputFiles.zip). 
An OSX-compatible executable and shell script are available [here](model_files/GMACS/tcsam02) and [here](model_files/GMACS/run_tcsam02_osx.sh). 
A Windows-compatible executable and batch file are available [here](model_files/GMACS/gmacs.exe) and [here](model_files/GMACS/run_tcsam02_bat.sh). 

To run the model: 

  * extract files in the zip file to a folder (the "top-level folder")
  * copy the gmacs executable (gmacs or gmacs.exe) to the top-level folder
    - you may have to set permissions to allow the executable to run
  * create a sub-folder "run1" and copy the shell script or batch file into it
    - you may have to set permissions to allow the shell script/batch file to run
  * run the shell script/batch file from a terminal/command window opened in the sub-folder
  
The ["gmacs.rep1"](model_files/GMACS/gmacs.rep1) is the main results file, although many other ancillary files are produced for a model run. 
If you have the necessary R packages installed, you can run the scripts in "r_CheckRep1.Run1a.R" and "r_CheckRep2.Run1a.R" to make plots comparing results 
from the bespoke model and GMACS. Note that ["rda_gmacs_rep.RData"](model_runs/GMACS/rda_gmacs_rep.RData) provides the version of "gmacs.rep1" already converted 
to R format.

## Hybrid Chionoecetes spp.

  * [Hybrid considerations](presentations/HybridConsiderations.pdf)

## Management Tier System and Proxies for $F_{MSY}$ and $B_{MSY}$

The NPFMC's management tier system for assessing Bering Sea and Aleutian Islands crab stocks is described [here](documents/2025.SAFE_Intro_TierSystem.pdf). The 
Tanner crab stock is assessed in Tier 3, which uses SPR-derived proxies for $F_{MSY}$ and $B_{MSY}$ to determine the federal overfishing limit (OFL) for the stock. 
For Tier 3, the proxies are $F_{MSY} = F_{35%}$ and $B_{MSY} = B_{35%}$. The [presentation](presentations/ReferencePointsAndTierConsiderations.pdf) will be used to 
motivate discussion on whether other SPR rates might be more appropriate for achieving sustainable harvest rates given a closer examination of Tanner crab life history.
 
## Spatial considerations and potential stock structure

  * [Spatial considerations](presentations/SpatialConsiderations.pdf)


