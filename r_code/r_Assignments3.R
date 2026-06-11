#--CIE Review assignments: Day 3
#--"Active Project" should be 2026-05_TannerCrab
#--working directory should be "this" file's folder
require(ggplot2);

#--load project info----
fn = file.path(rstudioapi::getActiveProject(),paste0("rda_project_info.",wtsUtilities::getOperatingSystem(),".RData"));
s = wtsUtilities::getObj(fn);

#--load 2025 assessment model----
lstT = wtsUtilities::getObj(s$results$Asmt2025)

##--extract model-estimated population biomass by x,m,z----
dfrPBfxmz = rTCSAM02::getMDFR.Pop.Biomass(lstT,cast="x+m+z") |> 
             dplyr::select(y,x,m,z,val) |> 
             dplyr::mutate(y=as.numeric(y),z=as.numeric(z),fleet="--",category="population",type="biomass");
##--extract model-estimated total catch mortality by x,m,z----
dfrCBfxmz = rTCSAM02::getMDFR.Fisheries.CatchBiomass(lstT,category="total mortality",cast="x+m+z") |> 
             dplyr::select(fleet,y,x,m,z,val) |> 
             dplyr::mutate(y=as.numeric(y),z=as.numeric(z),category="total catch mortality",type="biomass");
##--extract model-estimated retained catch mortality by f,x,m,z----
dfrRBfxmz = rTCSAM02::getMDFR.Fisheries.CatchBiomass(lstT,category="retained",cast="x+m+z") |> 
             dplyr::select(fleet,y,x,m,z,val) |> 
             dplyr::mutate(y=as.numeric(y),z=as.numeric(z),category="retained catch",type="biomass");
##--combine above
dfrPCBfxmz = dplyr::bind_rows(dfrPBfxmz,dfrCBfxmz,dfrRBfxmz);
readr::write_csv(dfrPCBfxmz,"Biomass_Pop-TotCatchMort-RetCatch.csv")
rm(dfrPBfxmz,dfrCBfxmz,dfrRBfxmz);

##--extract model-estimated population abundance by x,m,z----
dfrPAfxmz = rTCSAM02::getMDFR.Pop.Abundance(lstT,cast="x+m+z") |> 
             dplyr::select(y,x,m,z,val) |> 
             dplyr::mutate(y=as.numeric(y),z=as.numeric(z),fleet="--",category="population",type="abundance");
##--extract model-estimated total catch mortality by f,x,m,z----
dfrCAfxmz = rTCSAM02::getMDFR.Fisheries.CatchAbundance(lstT,category="total mortality",cast="x+m+z") |> 
             dplyr::select(fleet,y,x,m,z,val) |> 
             dplyr::mutate(y=as.numeric(y),z=as.numeric(z),category="total catch mortality",type="abundance");
##--extract model-estimated retained catch mortality by f,x,m,z----
dfrRAfxmz = rTCSAM02::getMDFR.Fisheries.CatchAbundance(lstT,category="retained",cast="x+m+z") |> 
             dplyr::select(fleet,y,x,m,z,val) |> 
             dplyr::mutate(y=as.numeric(y),z=as.numeric(z),category="retained catch",type="abundance");
##--combine above
dfrPCAxmz = dplyr::bind_rows(dfrPAfxmz,dfrCAfxmz,dfrRAfxmz);
rm(dfrPAfxmz,dfrCAfxmz,dfrRAfxmz);

##--extract model-estimated recruitment and spawning stock biomass----
scl=1;
dfrR_SSB = rTCSAM02::getMDFR.SdRep.RecAndSSB(lstT) |> 
             dplyr::select(variable,y,x,est,lci,uci) |> 
             dplyr::mutate(est=ifelse(variable=="lnRec",exp(est)/scl,est),
                           lci=ifelse(variable=="lnRec",exp(lci)/scl,lci),
                           uci=ifelse(variable=="lnRec",exp(uci)/scl,uci),
                           variable=ifelse(variable=="lnRec","R",paste(x,"MB")));
