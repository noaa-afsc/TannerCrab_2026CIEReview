# 2026 Tanner Crab CIE Review Repository
**NOTE: under construction**

This is the repository for materials for the 2026 CIE Review of the Bering Sea and Aleutian Islands (BSAI) Tanner crab assessment, which will be held June 9-11, 2026 at the Alaska Fisheries Science Center (Seattle, WA).
The BSAI Tanner crab (*Chionoecetes bairdi*) stock is one of [10 crab stocks](https://www.npfmc.org/fisheries/bsai-crab/) managed by the North Pacific Fishery Management Council ([NPFMC](https://www.npfmc.org/about-the-council/)). 
The review will include a number of topics including: 

  * assessment model transition from the current bespoke modeling framework (TCSAM02) to the framework used by the majority of 
BSAI crab stock assessments, the Generalized Model for Assessing Crustacean Stocks (GMACS)
  * issues related to whether, and if so how, hybrid *Chionoecetes* should be included in the assessment
  * issues related to the use of "side-by-side" survey data from gear calibration studies
  * issues related to management Tier level, management proxies for $F_{MSY}$ and $B_{MSY}$, and OFL and ABC calculations
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

### ADFG Fisheries



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

## TCSAM02 to GMACS bridging analysis

The Tanner crab assessment has been based on the Tanner Crab Stock Assessment Model, version
2, (TCSAM02) modeling framework for over 10 years. The framework
consists of an AD Model Builder C++ template, native C++ code in several additional source and include files, 
and the C++ library wtsADMB. Although the framework has many options for defining an implemented “model”, 
the characteristics of its current population dynamics are limited to *Chionoecetes* species (i.e., Tanner crab and snow crab) 
that exhibit a terminal molt to maturity 
and thus cannot be applied to, e.g., king crab species. Consequently, TCSAM02 does not provide a unified
modeling framework for Alaskan crab stocks in a manner similar to that of the Generalized Model for
Assessing Crustacean Stocks (GMACS). For several years, the SSC has requested that the Tanner
crab assessment be moved to the GMACS framework so that all BSAI crab stocks are assessed using
the framework, but it and the CPT have also indicated the necessity (CPT 2024, SSC 2024) that
a detailed bridging analysis from the bespoke model be conducted that demonstrates the resulting
GMACS model is “exactly equivalent” to the bespoke model before the transition from TCSAM02
to GMACS should occur for the Tanner crab assessment. Up to now, three of the major barriers
to development of an “exactly equivalent” GMACS Tanner crab model have been: 1) that several
features utilized in the bespoke model have not been available in GMACS, 2) that results were not
consistently reported between the two frameworks, and 3) that functions did not exist that allowed
comparison of results between the two models. From the September 2025 assessment to just prior to the May 2026 
CPT meeting, these barriers have, on the whole, been surmounted. In the [report](documents/TannerCrabReport_2026-05.pdf) (Section 4) 
submitted to the CPT in May, three areas were identified that remained to be better reconciled: 1) fitting to “extended” size composition data, 
2) calculation of the OFL and 3) priors and penalties applied to the model likelihood/objective function. 

Following the main report, the causes associated with the first issue, the differences between the two frameworks when fitting 
"extended" size compositions, were identified and GMACS code was updated with the result that extended size comps are now 
handled in identical fashion in both frameworks (for details, see the 
[addendum](documents/TannerCrabReport_2026-05_Addendum.pdf) to the May report to the CPT). Discussion at the May CPT meeting concluded 
that the ~8% difference in OFLs calculated by the two frameworks, the basis for the second issue,was probably not unreasonable given 
the different algorithms used to determine the OFL in the two frameworks. Thus, the remaining bridging issue is to reconcile the 
priors and penalties applied in GMACS with those applied in TCSAM02 for this stock and to demonstrate model convergence in GMACS "to" 
the solution found by the TCSAM02 model. Given that one cannot reasonably expect the two solutions to be identical in all respects, given 
finite numerical precision if nothing else, the question becomes: how close if "good enough" and what criteria should be used to 
determine "close".

## *Chionoecetes spp.* Hybrids

Abundance and biomass estimates for Chionoecetes spp. hybrids (“hybrids”) in the 2025 NMFS
EBS shelf bottom trawl survey were “unprecedented” across all sex/size/maturity categories (Zacher
et al. 2026). Male crab ≥ 78 mm CW had a biomass estimate of 37,068 t, representing a 471%
increase since 2024; it was also considerably higher than the previous 20-year average of 4,459 t.
Twenty percent of all the Chionoecetes spp. males ≥ 101 mm CW in the eastern Bering Sea were
hybrids. A peak of immature hybrid males was observed in the 50 – 90 mm range in the 2024 survey
and many of these crab appeared to reach mature and legal size classes for snow crab in 2025. In
addition, there were further large increases in juvenile size classes, portending future recruitment
to those size classes. Hybrid males were found in highest abundance on the middle shelf to the east
and northeast of the Pribilof Islands, with small male hybrids found further north than large male
hybrids.

The astonishing abundance of hybrids led to concerns expressed by stakeholders and reiterated by
the CPT and SSC at their respective Fall, 2025 meetings (CPT 2025, SSC 2025). Comments
from stakeholders included questions about the ways in which hybrids are included or excluded
from the Tanner crab OFL calculation, as well as the need for rapid action on developing a plan
for incorporating hybrids into assessments and regulations such that harvesters can take advantage
of the increasing abundance of hybrid males at industry-preferred sizes.

As a preliminary step, Tanner and snow crab stock authors were asked to complete three model
sensitivity runs (including size composition and abundance data):

  * With hybrids included in survey data;
  * With hybrids included in catch data;
  * With hybrids included in both survey and catch data

Results from this request are addressed in section 3 of the May 2026 Tanner crab [report](documents/TannerCrabReport_2026-05.pdf) 
(also referenced above) to the CPT. A link to the Tanner crab-specific portion of the presentation to the CPT is 
provided [here](presentations/HybridConsiderations.pdf). 
The entire presentation, including more information on survey results, fishery results, and results from snow crab model runs c
an be found [here](https://meetings.npfmc.org/CommentReview/DownloadFile?p=a1c76f9e-703a-477f-9128-447b1d6d6054.pdf&fileName=PPT_Hybrid%20topics_052026.pdf).

## Side-by-side selectivity studies

Side-by-side (SBS) trawling experiments were first conducted jointly by the National Marine
Fisheries’ Alaska Science Center (NMFS AFSC) and the Bering Sea Fisheries Research Foundation
(BSFRF) in the eastern Bering Sea (EBS) in 2009 and 2010 at standard NMFS EBS Shelf Survey stations in the northwest of the
EBS shelf survey grid to characterize the selectivity of the NMFS bottom trawl survey gear for
snow crab, *Chionoecetes opiliio* (Somerton et al. 2013). Similar studies were conducted in 2013,
2014, 2015, and 2016 in the southeastern EBS (Bristol Bay) to characterize the selectivity of the
NMFS bottom trawl survey gear for red king crab (*Paralithodes camtschaticus*) and Tanner crab
(*C. bairdi*), as well as further west on the EBS shelf in 2017 and 2018 to focus more
specifically on Tanner and snow crab.

In the current assessment model, annual survey biomass indices and size compositions from the BSFRF hauls are 
fit assuming they are absolute indices of abundance (i.e., *q* = 1) within the areas in which the BSFRF surveys were conducted. 
Incorporating these data into the assessment helps to determine the overall scale of the more comprehensive (in both space and time) 
relative indices provided by the NMFS EBS survey. An alternative approach to incorporating the BSFRF data into the assessment, 
rather than directly fitting the BSFRF biomass indices and size compositions in the model, has been to use the side-by-side, paired haul 
nature of the joint BSFRF-NMFS data to estimate the size selectivity of the NMFS gear outside the assessment model and use the 
results to inform survey selectivity within it. While not complete at this point, the CIE reviewers are asked to review the 
analysis to date and make recommendations on the steps needed to complete it ([link](documents/TannerCrab_SBS_Analysis.pdf)). Of note, 
the linked draft tech memo includes a haul-level analysis of the data for male Tanner crab but not for females; the 
latter follows an approach similar to that taken for males but the relevant section has not been written. 

A presentation on the side-by-side selectivity studies is available [here](presentations/SBS_Selectivity.pdf).

## Management Tier System and Proxies for $F_{MSY}$ and $B_{MSY}$

The NPFMC's management tier system for assessing Bering Sea and Aleutian Islands crab stocks is described [here](documents/2025.SAFE_Intro_TierSystem.pdf). The 
Tanner crab stock is assessed in Tier 3, which uses SPR-derived proxies for $F_{MSY}$ and $B_{MSY}$ to determine the federal overfishing limit (OFL) for the stock. 
For Tier 3, the proxies are $F_{MSY} = F_{35\\%}$ and $B_{MSY} = B_{35\\%}$. The [presentation](presentations/ReferencePointsAndTierConsiderations.pdf) will be used to 
motivate discussion on whether other SPR rates might be more appropriate for achieving sustainable harvest rates given a closer examination of Tanner crab life history.
 
## Spatial considerations and potential stock structure

Federal fisheries management treats Tanner crab in the EBS as a single stock: the stock assessment model 
fits survey and fishery data aggregated across the entire shelf and produces a single OFL and ABC for the entire area. 
In contrast, the State of Alaska sets separate TACs for the areas east and west of 166^o^W longitude 
based on area-specific harvest control rules. While Tanner crab in the EBS are regarded as a single stock, 
the population exhibits environmentally-driven and perhaps simply random annual changes on top of a fair degree of 
consistent or slowly-varying spatial structure together with annual changes. The reviewers are asked to consider the 
brief introduction presented in this [document](documents/TannerCrab_SpatialConsiderations.pdf) and consider whether future 
development of the assessment should quantitatively address spatial aspects of the stock. The figures are also provided in the 
[Spatial considerations](presentations/SpatialConsiderations.pdf) presentation.


# References

[CPT Report 2025](https://meetings.npfmc.org/CommentReview/DownloadFile?p=3506263a-4da3-4d77-aafb-ffbc694b69ef.pdf&fileName=C3%20CPT%20Report.pdf)

[SSC Report 2025](https://meetings.npfmc.org/CommentReview/DownloadFile?p=bb3958b1-c8ca-42ab-9391-c5673c6e9872.pdf&fileName=SSC%20Report%20Oct%202025_FINAL.pdf)

Somerton, D.A., Weinberg, K.L., and Goodman, S.E. 2013. Catchability of snow crab
(*Chionoecetes opilio*) by the eastern bering sea bottom trawl survey estimated using
a catch comparison experiment. Can. J. Fish. Aquat. Sci. 70: 1699–1708.
doi:dx.doi.org/10.1139/cjfas-2013-0100.

Thygesen,U.H., K. Kristensen, T. Jansen, J.E. Beyer. 2019. Intercalibration of survey methods using paired fishing operations and log-Gaussian Cox processes. 
ICES Journal of Marine Science, 76:1189–1199. https://doi.org/10.1093/icesjms/fsy191

Zacher, L.S., Hennessey, S.M., Richar, J.I., Fedewa, E.J., Ryznar, E.R., and Litzow, M.A. 2026. The
2025 eastern and northern Bering Sea continental shelf trawl surveys: Results for commercial
crab species. U.S. National Oceanic and Atmospheric Administration (NOAA). https://repository.library.noaa.gov/view/noaa/72630


