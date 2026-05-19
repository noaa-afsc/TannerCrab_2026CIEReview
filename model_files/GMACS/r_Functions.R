##--use wtsGMACS::readModelResults(folders) to read GMACS results into lstGMACS
#--extract ? for fits to indices
# rep = lstGs$repsLst$gmacs;
# dfr = rep$Index_fit_summary;

#--extract recruitment-----
extractRecruitment<-function(lstGMACS,lstTCSAM02=NULL){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractRecruitment(lstGMACS$repsLst[[nm]]) |> 
                    dplyr::mutate(case=nm);
    }
    dfrG = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    rep = lstGMACS;
    zBs = rep$size_midpoints;
    nZBs = length(zBs);
    #--get gmacs recruitment----
    dfrG = rep$R_y |> dplyr::select(y=year,val=est) |> 
             dplyr::mutate(case="gmacs",.before=1) |> 
             dplyr::mutate(y=as.numeric(y),
                           val=2*as.numeric(val));    #--gmacs recruitment is 0.5*tcsam recruitment
  }
  if (!is.null(lstTCSAM02)) {
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02)
    #--getTCSAM02 recruitment
    dfrT = rCompTCMs::extractMDFR.Pop.Recruitment(lst) |> 
             dplyr::select(case,y,val);
    dfrG = dplyr::bind_rows(dfrG,dfrT);
  }
  return(dfrG);
}
extractSizeAtRecruitment<-function(lstGMACS,lstTCSAM02=NULL){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractSizeAtRecruitment(lstGMACS$repsLst[[nm]]) |> 
                    dplyr::mutate(case=nm);
    }
    dfrG = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    rep = lstGMACS;
    #--get gmacs size at recruitment----
    dfrG = rep$R_z |> dplyr::select(x=sex,z=size,val=est) |> 
             dplyr::mutate(case="gmacs",.before=1) |> 
             dplyr::mutate(y="all",
                           z=as.numeric(z),
                           val=as.numeric(val));
  }
  if (!is.null(lstTCSAM02)) {
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02)
    #--getTCSAM02 recruitment
    dfrT = rCompTCMs::extractMDFR.Pop.RecSizeDistribution(lst) |> 
             dplyr::select(case,y,x,z,val);
    dfrG = dplyr::bind_rows(dfrG,dfrT);
  }
  return(dfrG);
}

#--extract cohort progression----
extractCohortProgression<-function(lstGMACS,lstTCSAM02=NULL,cast="x+m+s+z",gmacsType=2){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractCohortProgression(lstGMACS$repsLst[[nm]],NULL,cast=cast,gmacsType=gmacsType) |> 
                    dplyr::mutate(case=nm);
    }
    dfrG = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    rep = lstGMACS;
    zBs = rep$size_midpoints;
    nZBs = length(zBs);
    gcast = c("y",stringr::str_split_1(cast,"\\+"));
    gcast = rlang::syms(gcast);
    #--get gmacs recruitment----
    type = paste0("CohortProgression",gmacsType)
    dfrG = rep[[type]] |> 
             dplyr::filter(season=="1") |>
             tidyr::pivot_longer(cols=5+1:nZBs,names_to="z",values_to="val") |> 
             dplyr::select(y=year,x=sex,m=maturity,s=shell,z,val) |> 
             dplyr::mutate(y=as.numeric(y)-1,
                           s=paste(s,"shell"),
                           val=as.numeric(val)) |>
             dplyr::group_by(!!!gcast) |> 
             dplyr::summarize(val=sum(val)) |> 
             dplyr::ungroup() |> 
             dplyr::mutate(case="gmacs",.before=1); 
  }
  if (!is.null(lstTCSAM02)) {
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02)
    #--getTCSAM02 population abundance
    tcast = c("case","y",stringr::str_split_1(cast,"\\+"),"val")
    tcast = rlang::syms(tcast);
    dfrT = rCompTCMs::extractMDFR.Pop.CohortProgression(lst,cast=cast) |> 
             dplyr::select(!!!tcast);
    dfrG = dplyr::bind_rows(dfrG,dfrT);
  }
  return(dfrG)
}

#--extract population abundance----
extractPopAbd<-function(lstGMACS,lstTCSAM02=NULL,cast="x+m+s+z"){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractPopAbd(lstGMACS$repsLst[[nm]],cast=cast) |> 
                    dplyr::mutate(case=nm);
    }
    dfrG = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    rep = lstGMACS;
    zBs = rep$size_midpoints;
    nZBs = length(zBs);
    gcast = c("y",stringr::str_split_1(cast,"\\+"));
    gcast = rlang::syms(gcast);
    #--get gmacs recruitment----
    dfrG = rep$N_YXMSZ |> tidyr::pivot_longer(cols=4+1:nZBs,names_to="z",values_to="val") |> 
             dplyr::select(y=year,x=sex,m=maturity,s=shell,z,val) |> 
             dplyr::mutate(y=as.numeric(y),
                           s=paste(s,"shell"),
                           val=as.numeric(val)) |>
             dplyr::group_by(!!!gcast) |> 
             dplyr::summarize(val=sum(val)) |> 
             dplyr::ungroup() |> 
             dplyr::mutate(case="gmacs",.before=1); 
  }
  if (!is.null(lstTCSAM02)) {
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02)
    #--getTCSAM02 population abundance
    tcast = c("case","y",stringr::str_split_1(cast,"\\+"),"val")
    tcast = rlang::syms(tcast);
    dfrT = rCompTCMs::extractMDFR.Pop.Abundance(lst,cast=cast) |> 
             dplyr::select(!!!tcast);
    dfrG = dplyr::bind_rows(dfrG,dfrT);
  }
  return(dfrG)
}

#--extract population biomass----
extractPopBio<-function(lstGMACS,lstTCSAM02=NULL,cast="x+m+s+z"){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractPopBio(lstGMACS$repsLst[[nm]],cast=cast) |> 
                    dplyr::mutate(case=nm);
    }
    dfrG = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    rep = lstGMACS;
    zBs = rep$size_midpoints;
    nZBs = length(zBs);
    gcast = c("y",stringr::str_split_1(cast,"\\+"));
    gcast = rlang::syms(gcast);
    #--get gmacs population biomass----
    dfrG = rep$B_YXMSZ |> tidyr::pivot_longer(cols=4+1:nZBs,names_to="z",values_to="val") |> 
             dplyr::select(y=year,x=sex,m=maturity,s=shell,z,val) |> 
             dplyr::mutate(y=as.numeric(y),
                           s=paste(s,"shell"),
                           val=as.numeric(val)) |>
             dplyr::group_by(!!!gcast) |> 
             dplyr::summarize(val=sum(val)) |> 
             dplyr::ungroup() |> 
             dplyr::mutate(case="gmacs",.before=1); 
  }
  if (!is.null(lstTCSAM02)) {
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02)
    #--getTCSAM02 recruitment
    tcast = c("case","y",stringr::str_split_1(cast,"\\+"),"val")
    tcast = rlang::syms(tcast);
    dfrT = rCompTCMs::extractMDFR.Pop.Biomass(lst,cast=cast) |> 
             dplyr::select(!!!tcast);
    dfrG = dplyr::bind_rows(dfrG,dfrT);
  }
  return(dfrG)
}

#--extract mean growth rates
extractMeanGrowth<-function(lstGMACS,lstTCSAM02=NULL){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractMeanGrowth(lstGMACS$repsLst[[nm]]) |> 
                    dplyr::mutate(case=nm);
    }
    dfr = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    rep = lstGMACS;
    #--get gmacs mean post-molt size by premolt size----
    dfr = rep$`Mean growth` |> 
            dplyr::select(x=sex,z=premolt_size,val=mean_postmolt_size) |>
                  dplyr::mutate(case=paste("gmacs"),
                                z=as.numeric(z),
                                val=as.numeric(val)) |>
                  rCompTCMs::getMDFR.CanonicalFormat();
  }
  if (!is.null(lstTCSAM02)){
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02)
    dfrT = rCompTCMs::extractMDFR.Pop.MeanGrowth(lst);
    dfr = dplyr::bind_rows(dfr,dfrT);
  }
  return(dfr);
}

#--extract growth matrices
extractGrowthMatrices<-function(lstGMACS,lstTCSAM02=NULL){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractGrowthMatrices(lstGMACS$repsLst[[nm]]) |> 
                    dplyr::mutate(case=nm);
    }
    dfr = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    rep = lstGMACS;
    #--get gmacs growth matrices----
    dfr = rep$`growth_matrix` |> 
            tidyr::pivot_longer(4:ncol(rep$`growth_matrix`)) |> 
            dplyr::select(x=sex,pc=block,z=premolt_size,zp=name,val=value) |>
                  dplyr::mutate(case=paste("gmacs"),
                                m="immature",
                                z=as.numeric(z),
                                zp=as.numeric(zp),
                                val=as.numeric(val)) |>
                  rCompTCMs::getMDFR.CanonicalFormat();
  }
  if (!is.null(lstTCSAM02)){
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02)
    dfrT = rCompTCMs::extractMDFR.Pop.GrowthMatrices(lst);
    dfr = dplyr::bind_rows(dfr,dfrT);
  }
  return(dfr);
}

#--extract probability of molt to maturity
extractPrM2M<-function(lstGMACS,lstTCSAM02=NULL){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractPrM2M(lstGMACS$repsLst[[nm]]) |> 
                    dplyr::mutate(case=nm);
    }
    dfr = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    rep = lstGMACS;
    #--get maturation probabilities----
    dfr = rep$`prMature` |> 
            tidyr::pivot_longer(3:ncol(rep$`prMature`)) |> 
            dplyr::select(x=sex,y=year,z=name,val=value) |>
                  dplyr::mutate(case=paste("gmacs"),
                                m="immature",
                                z=as.numeric(z),
                                val=as.numeric(val)) |>
                  rCompTCMs::getMDFR.CanonicalFormat();
  }
  if (!is.null(lstTCSAM02)){
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02)
    dfrT = rCompTCMs::extractMDFR.Pop.PrM2M(lst);
    if ((class(dfrT$y)=="character")||(class(dfr$y)=="character")){
      dfrT$y = as.character(dfrT$y);
      dfr$y  = as.character(dfr$y);
    }
    dfr = dplyr::bind_rows(dfr,dfrT);
  }
  return(dfr);
}

#--extract fishing mortality rates----
extractFisheryFs<-function(lstGMACS,lstTCSAM02=NULL,season=2){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractFisheryFs(lstGMACS$repsLst[[nm]],season=season);
      if (!is.null(lst[[nm]])) lst[[nm]] = lst[[nm]] |> dplyr::mutate(case=nm);
    }
    dfr = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    dfrp = lstGMACS$`Fully-selected_FM_by_season_sex_and_fishery`;
    if (is.null(dfrp)) {
      dfrp = lstGMACS$`Fully-selected_capture_rate_by_season_sex_and_fishery`;
      if (is.null(dfrp)) return(NULL);
    }
    #--get gmacs fishery Fs----
    dfr = dfrp |> # chkNames() |>
            dplyr::select(c(1:3,3+season));
    dfr$val = as.numeric(dfr[[as.character(season)]]);
    dfr  = dfr |> dplyr::select(fleet,y=year,x=sex,val) |>
                  dplyr::mutate(y=as.numeric(y),
                                case=paste("gmacs")) |>
                  rCompTCMs::getMDFR.CanonicalFormat();
  }
  if (!is.null(lstTCSAM02)){
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02)
    dfrT = rCompTCMs::extractMDFR.Fisheries.Catchability(lst);
    dfr = dplyr::bind_rows(dfr,dfrT);
  }
  return(dfr);
}

#--extract total fishing mortality----
extractTotalFishingMortality<-function(lstGMACS,lstTCSAM02=NULL){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractTotalFishingMortality(lstGMACS$repsLst[[nm]]) |> 
                    dplyr::mutate(case=nm);
    }
    dfr = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    rep = lstGMACS;
    dfr = rep$Summary |> dplyr::select(y=year,val=tot_mortality) |> 
             dplyr::mutate(y=as.numeric(y),
                           val=as.numeric(val),
                           case="gmacs");
  }
  if (!is.null(lstTCSAM02)){
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02)
    dfrT = rCompTCMs::extractMDFR.Fisheries.CatchBiomass(lst,
                                                         category="total mortality",
                                                         cast="y") |> 
             dplyr::filter(type=="predicted") |> 
             dplyr::group_by(case,y) |> 
             dplyr::summarize(val=sum(val,na.rm=TRUE)) |> 
             dplyr::ungroup() |>
             dplyr::select(case,y,val);
    dfr = dplyr::bind_rows(dfr,dfrT); 
  }
  return(dfr)
}


#--extract retained catch mortality----
extractRetainedCatchMortality<-function(lstGMACS,lstTCSAM02=NULL){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractRetainedCatchMortality(lstGMACS$repsLst[[nm]]) |> 
                    dplyr::mutate(case=nm);
    }
    dfr = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    rep = lstGMACS;
    dfr = rep$Summary |> dplyr::select(y=year,val=ret_mortality) |> 
             dplyr::mutate(y=as.numeric(y),
                           val=as.numeric(val),
                           case="gmacs");
  }
  if (!is.null(lstTCSAM02)){
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02)
      dfrT = rCompTCMs::extractMDFR.Fisheries.CatchBiomass(list(tcsam=lstT),
                                                           category="retained",
                                                           cast="y") |> 
               dplyr::filter(type=="predicted",category=="retained") |> 
               dplyr::select(case,y,val);
      dfr = dplyr::bind_rows(dfr,dfrT); 
  }
  return(dfr)
}

#--extract discard catch mortality----
extractDiscardCatchMortality<-function(lstGMACS,lstTCSAM02=NULL){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractDiscardCatchMortality(lstGMACS$repsLst[[nm]]) |> 
                    dplyr::mutate(case=nm);
    }
    dfr = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    rep = lstGMACS;
    flts = rep$Fleets;
    tdys = c("year",flts,"ret_mortality");
    # dfr = rep$Summary |> dplyr::select(y=year,TCF,SCF,RKF,GFA,ret_mortality) |> 
    #          dplyr::mutate(dplyr::across(c(y,TCF,SCF,RKF,GFA),as.numeric)) |> 
    #          dplyr::mutate(TCF=as.numeric(TCF)-as.numeric(ret_mortality)) |> dplyr::select(!ret_mortality) |> 
    #          tidyr::pivot_longer(c(TCF,SCF,RKF,GFA),names_to="fleet",values_to="val") |> 
    #          dplyr::mutate(case="gmacs");
    dfr = rep$Summary |> dplyr::select(tidyselect::all_of(tdys)) |> 
             dplyr::mutate(dplyr::across(tidyselect::all_of(tdys),as.numeric)) |> 
             dplyr::mutate(y=year,
                           TCF=TCF-ret_mortality) |> 
             dplyr::select(!c(ret_mortality,year)) |> 
             tidyr::pivot_longer(tidyselect::all_of(flts),names_to="fleet",values_to="val") |> 
             dplyr::mutate(case="gmacs");
  }
  if (!is.null(lstTCSAM02)){
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02)
      dfrT = rCompTCMs::extractMDFR.Fisheries.CatchBiomass(list(tcsam=lstT),
                                                           category="discard mortality",
                                                           cast="y") |> 
               dplyr::filter(type=="predicted",category=="discard mortality") |> 
               dplyr::select(case,fleet,y,val);
      dfr = dplyr::bind_rows(dfr,dfrT); 
  }
  return(dfr)
}

#--extract fishery capture biomass----
extractFisheryCaptureBiomass<-function(lstGMACS,lstTCSAM02=NULL){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractFisheryCaptureBiomass(lstGMACS$repsLst[[nm]]);
      if (!is.null(lst[[nm]])) lst[[nm]] = lst[[nm]] |>  dplyr::mutate(case=nm);
    }
    dfr = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    dfrp = lstGMACS$`Predicted_capture_biomass-at-size`;
    if (is.null(dfrp)) return(NULL);
    dfr = dfrp |>
             tidyr::pivot_longer(4:ncol(dfrp),
                                 names_to="z",values_to="val") |> 
             dplyr::rename(y=year,x=sex) |> 
             dplyr::mutate(y=as.numeric(y),
                           z=as.numeric(z),
                           val=as.numeric(val)) |> 
             dplyr::group_by(y,fleet,x) |> 
             dplyr::summarize(val=sum(val,na.rm=TRUE)) |> 
             dplyr::ungroup() |> 
             dplyr::mutate(case="gmacs");
  }
  if (!is.null(lstTCSAM02)){
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02);
      dfrT = rCompTCMs::extractMDFR.Fisheries.CatchBiomass(list(tcsam=lstT),
                                                           category="captured",
                                                           cast="y+x") |> 
               dplyr::filter(type=="predicted",category=="captured") |> 
               dplyr::select(case,fleet,y,x,val);
      dfr = dplyr::bind_rows(dfr,dfrT); 
  }
  return(dfr)
}

#--extract fishery capture abundance----
extractFisheryCaptureAbundance<-function(lstGMACS,lstTCSAM02=NULL){
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractFisheryCaptureAbundance(lstGMACS$repsLst[[nm]]);
      if (!is.null(lst[[nm]])) lst[[nm]] = lst[[nm]] |>  dplyr::mutate(case=nm);
    }
    dfr = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    dfrp = lstGMACS$`Predicted_capture_abundance-at-size`;
    if (is.null(dfrp)) return(NULL);
    dfr = dfrp |>
             tidyr::pivot_longer(4:ncol(dfrp),
                                 names_to="z",values_to="val") |> 
             dplyr::rename(y=year,x=sex) |> 
             dplyr::mutate(y=as.numeric(y),
                           z=as.numeric(z),
                           val=as.numeric(val)) |> 
             dplyr::group_by(y,fleet,x) |> 
             dplyr::summarize(val=sum(val,na.rm=TRUE)) |> 
             dplyr::ungroup() |> 
             dplyr::mutate(case="gmacs");
  }
  if (!is.null(lstTCSAM02)){
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02);
      dfrT = rCompTCMs::extractMDFR.Fisheries.CatchAbundance(list(tcsam=lstT),
                                                           category="captured",
                                                           cast="y+x") |> 
               dplyr::filter(type=="predicted",category=="captured") |> 
               dplyr::select(case,fleet,y,x,val);
      dfr = dplyr::bind_rows(dfr,dfrT); 
  }
  return(dfr)
}

#--extract predicted survey biomass----
extractPredictedSurveyBiomass<-function(lstGMACS,lstTCSAM02=NULL,fleetG=NULL,fleetT=fleetG){
  if (is.null(fleetG)) stop("must supply a survey fleet name");
  if (inherits(lstGMACS,"gmacs_reslst")) {
    lst=list();
    for (nm in names(lstGMACS$repsLst)){
      #--testing: nm = names(lstGMACS$repsLst)[1];
      lst[[nm]] = extractPredictedSurveyBiomass(lstGMACS$repsLst[[nm]],fleetG=fleetG);
      if (!is.null(lst[[nm]])) lst[[nm]] = lst[[nm]] |>  dplyr::mutate(case=nm);
    }
    dfr = dplyr::bind_rows(lst);
    rm(lst);
  } else if (inherits(lstGMACS,"gmacs_rep1")) {
    dfrp = lstGMACS$`Index_fit_summary`;
    if (is.null(dfrp)) return(NULL);
    dfr = dfrp |>
           dplyr::filter(stringr::str_starts(fleet,fleetG),units %in% c("biomass","biommass")) |> 
           dplyr::select(y=year,x=sex,m=maturity,val=prd) |> 
           dplyr::mutate(y=as.numeric(y),val=as.numeric(val),m=ifelse(m=="undetermined","all",m)) |> 
           dplyr::mutate(case="gmacs");
  }
  if (!is.null(lstTCSAM02)){
    if (inherits(lstTCSAM02,"tcsam02.resLst")) lst = list(tcsam=lstTCSAM02);
      dfrT1 = rCompTCMs::extractMDFR.Surveys.Biomass(lst,cast="x+m") |> 
               dplyr::filter(stringr::str_starts(fleet,fleetT),x=="female") |> 
               dplyr::select(case,y,x,m,val);
      dfrT2 = rCompTCMs::extractMDFR.Surveys.Biomass(lst,cast="x") |> 
               dplyr::filter(stringr::str_starts(fleet,fleetT),x=="male") |> 
               dplyr::select(case,y,x,m,val);
      dfr = dplyr::bind_rows(dfr,dfrT1,dfrT2); 
  }
  return(dfr)
}

# #--extract survey Qs (INCOMPLETE)--DON'T USE----
# extractSurveyQs<-function(resGMACS,resTCSAM02=NULL){
#   if (inherits(lstGMACS,"gmacs_rep1")) rep = lstGMACS;
#   zBs = rep$size_midpoints;
#   nZBs = length(zBs);
#   #--get gmacs survey Qs----
#   dfrM = rep$`Natural_mortality-by-class` |> chkNames() |>
#                  tidyr::pivot_longer(cols=4+(1:nZBs),names_to="z",values_to="val") |>
#                  dplyr::rename(y=year,x=sex,m=maturity,s=shell_cond) |>
#                  dplyr::mutate(dplyr::across(c(1,5,6),as.numeric),
#                                case=paste("gmacs",names(resGMACS)[1])) |>
#                  rCompTCMs::getMDFR.CanonicalFormat();
#   if (!is.null(resTCSAM02)){
#     if (inherits(lstTCSAM02,"tcsam02.resLst")) res = ;
#     dfrMt = rCompTCMs::extractMDFR.Pop.NaturalMortality(resTCSAM02,type="M_yxmsz") |>
#               dplyr::mutate(z=as.numeric(z));
#     if (length(unique(dfrM$x))==1) {
#       cols = names(dfrMt)[!(names(dfrMt) %in% c("x","val"))];
#       dfrMt = dfrMt |> dplyr::mutate(x="undetermined") |>
#                 dplyr::group_by(!!!rlang::syms(cols)) |>
#                 dplyr::summarize(val=mean(val)) |>
#                 dplyr::ungroup() |>
#                 dplyr::mutate(x="undetermined");
#     }
#     if (length(unique(dfrM$m))==1) {
#       cols = names(dfrMt)[!(names(dfrMt) %in% c("m","val"))];
#       dfrMt = dfrMt |> dplyr::mutate(m="undetermined") |>
#                 dplyr::group_by(!!!rlang::syms(cols)) |>
#                 dplyr::summarize(val=mean(val)) |>
#                 dplyr::ungroup() |>
#                 dplyr::mutate(m="undetermined");
#     }
#     if (length(unique(dfrM$s))==1) {
#       cols = names(dfrMt)[!(names(dfrMt) %in% c("s","val"))];
#       dfrMt = dfrMt |> dplyr::mutate(s="undetermined") |>
#                 dplyr::group_by(!!!rlang::syms(cols)) |>
#                 dplyr::summarize(val=mean(val)) |>
#                 dplyr::ungroup() |>
#                 dplyr::mutate(s="undetermined");
#     }
#     dfrM = dplyr::bind_rows(dfrM,dfrMt);
#   }
#   return(dfrM);
# }

# #--compare selectivity Qs (INCOMPLETE: DON'T USE)----
# compareSurveyQs<-function(dfr,ylab="Survey Q",plotZ=FALSE){
#   ##--drop unnecessary columns----
#   cols = names(dfr)[names(dfr) %in% c("case","y","x","m","s","z","val")];
#   ncls = length(cols);
#   tmp0 = dfr |> dplyr::select(tidyselect::all_of(cols));
#   ##--collect values and create labels
#   tmp1 = wtsUtilities::collectValuesByGroup(tmp0,collect="y",names_from="z",values_from="val");
#   tmp1 = tmp1 |> dplyr::rowwise() |>
#                  dplyr::mutate(ylab=wtsUtilities::collapseIntegersToString(y),.before=y) |>
#                  dplyr::ungroup();
#   ##--cross cases----
#   tmp2  = tmp1 |> dplyr::select(tidyselect::any_of(1:ncls));
#   tmp2a = tmp2 |> dplyr::cross_join(tmp2);
#   tmp3  = tmp2a |> dplyr::filter(!(case.x==case.y),(stringr::str_starts(case.x,"gmacs")));
#   if (nrow(tmp3)>0){
#     ##--select unique case combinations out of duplicates and self-crosses
#     tmp4 = tmp3 |> dplyr::rowwise() |> 
#              dplyr::mutate(check=(any(unlist(y.y) %in% unlist(y.x)))) |> 
#              dplyr::ungroup() |> dplyr::filter(check) |> 
#              dplyr::mutate(grp2=paste(group.x,group.y),
#                            ylabs=paste0(ylab.x,"\n",ylab.y),
#                            .before=1);
#     xcols = names(tmp4)[stringr::str_ends(names(tmp4),".x")]
#     ycols = names(tmp4)[stringr::str_ends(names(tmp4),".y")]
#     tmp5g = tmp4 |> dplyr::select(tidyselect::any_of(c("grp2","ylabs",xcols)));
#     names(tmp5g) = c("grp2","ylabs",stringr::str_remove(xcols,".x"));
#     tmp5t = tmp4 |> dplyr::select(tidyselect::any_of(c("grp2","ylabs",ycols)));
#     names(tmp5t) = c("grp2","ylabs",stringr::str_remove(ycols,".y"));
#     tmp5 = dplyr::bind_rows(tmp5g,tmp5t);
#   } else {
#     paste("got here");
#     tmp5 = tmp2 |> dplyr::mutate(grp2=group,
#                                   ylabs=ylab,
#                                   .before=1);
#   }
#   ##--convert to long format, with numeric z's, drop "y" column (don't need it)
#   ncls5 = ncol(tmp5);
#   tmp6 = tmp5 |> dplyr::inner_join(tmp1) |>
#            tidyr::pivot_longer(ncls5+1:(ncol(tmp1)-ncls),names_to="z",values_to="val") |>
#            dplyr::mutate(z=as.numeric(z)) |>
#            dplyr::select(!y);
#   
#   ##--make plot----
#   dfrp = tmp6 |> dplyr::distinct(grp2,ylabs,ylab,case,x,m,s,val);
#   if (length(unique(dfrp$s)) > 1) {
#     dfrp = dfrp |> dplyr::mutate(xms=paste0(m," ",x,"\n",s));
#   } else {
#     dfrp = dfrp |> dplyr::mutate(xms=paste0(m," ",x));
#   }
#   mx = max(max(dfrp$val),0.23);
#   p1 = ggplot(dfrp |> dplyr::filter(m=="immature"),
#               aes(x=xms,y=val,colour=case,fill=case)) + 
#          geom_col(position="dodge") + 
#          geom_hline(yintercept=0.23,linetype=3) + 
#          geom_hline(yintercept=0.00,linetype=3) + 
#          scale_y_continuous(limits=c(0,mx)) + 
#          facet_wrap(~ylab,ncol=1) + 
#          labs(x="Population class",y=ylab) + 
#          wtsPlots::getStdTheme() + 
#          theme(axis.title.x=element_blank(),
#                legend.direction="horizontal",
#                legend.position="inside",
#                legend.position.inside=c(0.01,0.99),
#                legend.justification.inside=c(0,1));
#   p2 = p1 %+% (dfrp |> dplyr::filter(m=="mature")) + theme(legend.position="none");
#   pg = cowplot::plot_grid(p1,p2,ncol=1,rel_heights=c(1,1.8))
#   return(pg)
# }
