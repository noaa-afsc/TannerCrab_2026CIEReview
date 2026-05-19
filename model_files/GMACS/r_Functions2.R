####--fits to aggregated index (survey) time series data----
getFitsToIndexTimeSeries<-function(lstGs,lstT=NULL,type="biomass"){
  lstGSDs = list();
  for (case_ in names(lstGs$repsLst)){
    dfrIFS = lstGs$repsLst[[case_]]$Index_fit_summary;
    lstGSDs[[case_]] = dfrIFS |> dplyr::filter(units==type) |> 
                        dplyr::select(y=year,fleet,x=sex,m=maturity,
                                      observed=obs,cv=actual_CV,predicted=prd,zscore=prsn_res) |> 
                        dplyr::mutate(dplyr::across(c(y,observed,cv,predicted,zscore),as.numeric),
                                      sdobs=sqrt(log(1+cv^2)),
                                      nll=0.5*log(2*pi*sdobs^2) + (0.5*zscore^2),
                                      residual=log(observed)-log(predicted)) |> 
                        dplyr::mutate(case=case_,.before=1);
  }
  if (!is.null(lstT)){ 
    if (type=="biomass"){
      dfrT = rTCSAM02::getMDFR.AllScores.Biomass(lstT,
                                                 fleet.type="survey",
                                                 catch.type="index");
    } else  {
      dfrT = rTCSAM02::getMDFR.AllScores.Abundance(lstT,
                                                 fleet.type="survey",
                                                 catch.type="index");
    }
    dfrT = dfrT |> 
            tidyr::pivot_wider(names_from=type,values_from=val) |> 
            dplyr::filter(observed > 0, !stringr::str_starts(.data$fleet,"SBS NMFS")) |> 
            dplyr::mutate(cv=sqrt(exp(sdobs^2)-1),
                          residual=log(observed)-log(predicted),
                          zscore=`z-score`,
                          nll=0.5*log(2*pi*sdobs^2) + (0.5*zscore^2),
                          m=ifelse(stringr::str_starts(m,"all"),"undetermined",m),
                          fleet=ifelse(stringr::str_starts(fleet,"NMFS"),"NMFS",fleet),
                          fleet=ifelse(stringr::str_starts(fleet,"SBS BSFRF"),"BSFRF",fleet),
                          case="tcsam") |> 
            dplyr::select(case,fleet,y,x,m,observed,cv,predicted,residual,zscore,nll);
    lstGSDs[["tscam"]] = dfrT;
  }
  dfr = dplyr::bind_rows(lstGSDs) |> 
          dplyr::mutate(xm=paste(m,x),
                        lwr=qlnorm(0.05,log(observed),sdlog=sdobs),
                        upr=qlnorm(0.95,log(observed),sdlog=sdobs));
  return(dfr);
}
plotFitsToIndexTSs <- function(dfr,f_="NMFS"){
                          dfrp = dfr |> dplyr::filter(fleet==f_);
                          p = ggplot(dfrp,aes(x=y,y=predicted,colour=case)) + 
                                geom_ribbon(aes(ymin=lwr,ymax=upr),data=dfrp |> dplyr::filter(!is.na(lwr)),fill="light gray",colour=NA) + 
                                geom_point(aes(y=observed)) + 
                                geom_line() + 
                                geom_hline(yintercept=0,linetype=3) + 
                                facet_grid(xm~.,scales="free_y") + 
                                wtsPlots::getStdTheme() + wtsPlots::noXT() + 
                                theme(legend.position="inside",
                                      legend.position.inside=c(1,1),
                                      legend.justification=c(1,1));
                          return(p);
                        }
# dfr = getFitsToIndexTimeSeries(lstGs,lstT) 
# p1 = plotFitsToIndexTSs(dfr,f="NMFS");      print(p1);
# p2 = p1 + scale_y_log10();                  print(p2);


####--fits to aggregated size comps----
getFitsToZCs<-function(lstGs,lstT=NULL){
  lstGSDs = list();
  for (case_ in names(lstGs$repsLst)){
    dfrIFS = lstGs$repsLst[[case_]]$Size_fit_summary;
    lstGSDs[[case_]] = dfrIFS |> 
                        dplyr::select(comp_type,y=year,fleet,x=sex,m=maturity,s=shell,z=size,
                                      inpSS=inpSS,observed=aggObs,predicted=aggPrd,residual=aggRes) |> 
                        dplyr::mutate(dplyr::across(c(y,z,inpSS,observed,predicted,residual),as.numeric)) |> 
                        dplyr::mutate(case=case_,.before=1);
  }
  if (!is.null(lstT)){
    lstT1 = list(tcsam=lstT);
    dfrT = dplyr::bind_rows(
               rCompTCMs::extractFits.SizeComps(lstT1,fleet.type="survey", catch.type="index")    |> dplyr::mutate(comp_type="total"),
               rCompTCMs::extractFits.SizeComps(lstT1,fleet.type="fishery",catch.type="retained") |> dplyr::mutate(comp_type="retained"),
               rCompTCMs::extractFits.SizeComps(lstT1,fleet.type="fishery",catch.type="total")    |> dplyr::mutate(comp_type="total"),
               rCompTCMs::extractFits.ZScores.PrNatZ(lstT1,fleet.type="survey", catch.type="index",   residuals.type="pearsons") |> dplyr::mutate(comp_type="total"),
               rCompTCMs::extractFits.ZScores.PrNatZ(lstT1,fleet.type="fishery",catch.type="retained",residuals.type="pearsons") |> dplyr::mutate(comp_type="retained"),
               rCompTCMs::extractFits.ZScores.PrNatZ(lstT1,fleet.type="fishery",catch.type="total",   residuals.type="pearsons") |> dplyr::mutate(comp_type="total")
            ) |> dplyr::mutate(val=ifelse((sign==">0")|(is.na(sign)),val,-1*val)) |> 
            dplyr::select(case,comp_type,y,fleet,x,m,s,z,type,val) |> 
            tidyr::pivot_wider(names_from="type",values_from="val") |> 
            dplyr::filter(!stringr::str_starts(.data$fleet,"SBS NMFS")) |> 
            dplyr::mutate(m=ifelse(stringr::str_starts(.data$m,"all"),"undetermined",m),
                          s=ifelse(stringr::str_starts(.data$s,"all"),"undetermined",s),
                          fleet=ifelse(stringr::str_starts(.data$fleet,"NMFS"),"NMFS",fleet),
                          fleet=ifelse(stringr::str_starts(.data$fleet,"SBS "),"BSFRF",fleet),
                          residual=pearsons);
    lstGSDs[["tscam"]] = dfrT;
  }
  dfr = dplyr::bind_rows(lstGSDs) |> 
          dplyr::mutate(fleet=stringr::str_replace_all(.data$fleet,"_"," "),
                        xm=paste(m,x));
  return(dfr);
}

plotFitsToZCs<-function(dfr){
  p = ggplot(dfr,aes(x=z,y=predicted,color=case)) + 
        geom_point(aes(y=observed),data=dfr,size=1) +  #--data
        geom_point(data=dfr |> dplyr::filter(case=="g1"),size=1) +        #--g1 pred as pts
        geom_line() +                                                     #--pred as line
        facet_wrap(~y) + 
        labs(x="size (mm CW)") + 
        wtsPlots::getStdTheme();
  return(p);
}

# dfr = getFitsToZCs(lstGs,lstT);
# plotFitsToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x=="male",  m=="undetermined",s=="undetermined"))
# plotFitsToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x=="female",m=="mature",      s=="undetermined"))
# plotFitsToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x=="female",m=="immature",    s=="undetermined"))
# plotFitsToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="retained",x=="male",  m=="undetermined",s=="undetermined"))
# plotFitsToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
# plotFitsToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
# plotFitsToZCs(dfr |> dplyr::filter(fleet=="SCF", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
# plotFitsToZCs(dfr |> dplyr::filter(fleet=="SCF", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
# plotFitsToZCs(dfr |> dplyr::filter(fleet=="RKF", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
# plotFitsToZCs(dfr |> dplyr::filter(fleet=="RKF", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
# plotFitsToZCs(dfr |> dplyr::filter(fleet=="GF All", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
# plotFitsToZCs(dfr |> dplyr::filter(fleet=="GF All", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))

plotDataToZCs<-function(dfr){
  p = ggplot(dfr,aes(x=z,y=observed,color=case)) + 
        geom_point(data=dfr |> dplyr::filter(case=="g1"),size=1) +        #--g1 pred as pts
        geom_line() +                                                     #--pred as line
        facet_wrap(~y) + 
        labs(x="size (mm CW)") + 
        wtsPlots::getStdTheme();
  return(p);
}
# plotDataToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x==  "male",m=="undetermined",s=="undetermined"))
# plotDataToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x=="female",m=="immature",    s=="undetermined"))
# plotDataToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x=="female",m==  "mature",    s=="undetermined"))
# plotDataToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="total",   x==  "male",m=="undetermined",s=="undetermined"))
# plotDataToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
# dfrp = dfr |> dplyr::filter(fleet=="TCF", comp_type=="total") |> 
#        dplyr::group_by(y,case) |> dplyr::summarize(observed=sum(observed),predicted=sum(predicted)) |> dplyr::ungroup();
# 
# dfrp = dfr |> dplyr::filter(fleet=="TCF", comp_type=="retained",x=="male",  m=="undetermined",s=="undetermined", y==2018);
# View(dfrp |> dplyr::select(case,z,observed) |> tidyr::pivot_wider(names_from="case",values_from="observed"))
# View(dfrp |> dplyr::select(case,z,predicted) |> tidyr::pivot_wider(names_from="case",values_from="predicted"))
# 
# dfrp = dfr |> dplyr::filter(fleet=="TCF", comp_type=="total",y==2022) |> 
#         dplyr::group_by(case,x) |> dplyr::summarize(observed=sum(observed),predicted=sum(predicted)) |> dplyr::ungroup();
# 
# dfrIFS |> dplyr::filter(fleet=="TCF",comp_type=="total",as.numeric(year) %in% 1991:1996,as.numeric(size)==127) |> 
#           dplyr::arrange(year) |> dplyr::group_by(year) |> dplyr::summarize(inpSS=sum(as.numeric(inpSS))) |> dplyr::ungroup();

plotResidualsToZCs<-function(dfr){
  p = ggplot(dfr,aes(x=z,y=residual,color=case,shape=case)) + 
        geom_line() + geom_point(size=1) + 
        facet_wrap(~y) + 
        labs(x="size (mm CW)") + 
        wtsPlots::getStdTheme();
  return(p);
}
# plotResidualsToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x==  "male",m=="undetermined",s=="undetermined"))
# plotResidualsToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x=="female",m=="mature",      s=="undetermined"))
# plotResidualsToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x=="female",m=="immature",    s=="undetermined"))
# plotResidualsToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="retained",x=="male",  m=="undetermined",s=="undetermined"))
# plotResidualsToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
# plotResidualsToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
# plotResidualsToZCs(dfr |> dplyr::filter(fleet=="SCF", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
# plotResidualsToZCs(dfr |> dplyr::filter(fleet=="SCF", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
# plotResidualsToZCs(dfr |> dplyr::filter(fleet=="RKF", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
# plotResidualsToZCs(dfr |> dplyr::filter(fleet=="RKF", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
# plotResidualsToZCs(dfr |> dplyr::filter(fleet=="GF_All", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
# plotResidualsToZCs(dfr |> dplyr::filter(fleet=="GF_All", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
# 
# dfre = dfr |> dplyr::select(case,fleet,comp_type,y,x,z,observed,predicted) |> 
#          dplyr::filter(fleet=="TCF",comp_type=="total") |> dplyr::select(case,y,x,z,observed,predicted) |>
#          tidyr::pivot_longer(c(observed,predicted),names_to="type",values_to="val") |> 
#          tidyr::pivot_wider(names_from=c("x","z"),values_from="val",names_sep="_") |> 
#          tidyr::pivot_longer(cols=4:67,names_to="xz",values_to="val") |> 
#          dplyr::mutate(xz=factor(xz,levels=unique(.data$xz)));
# 
# p = ggplot(mapping=aes(x=as.numeric(xz),y=val,color=case)) + 
#       geom_line() +facet_wrap(~y) +
#       wtsPlots::getStdTheme() + theme(axis.text.x=element_blank());
# p + (dfre |> dplyr::filter(type=="observed"))         
# p + (dfre |> dplyr::filter(type=="predicted"))         
# 
# p = ggplot(mapping=aes(x=as.numeric(xz),y=val,color=type)) + 
#       geom_line() +facet_wrap(~y) +
#       wtsPlots::getStdTheme() + theme(axis.text.x=element_blank());
# p + (dfre |> dplyr::filter(case=="tcsam"))         
# p + (dfre |> dplyr::filter(case=="g1"))         

  
  
  
  