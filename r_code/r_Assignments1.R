#--CIE Review assignments: Day 1
#--"Active Project" should be 2026-05_TannerCrab
require(ggplot2);

#--load project info----
fn = file.path(rstudioapi::getActiveProject(),paste0("rda_project_info.",wtsUtilities::getOperatingSystem(),".RData"));
s = wtsUtilities::getObj(fn);

#--load 2025 assessment model----
lstT = wtsUtilities::getObj(s$results$Asmt2025)

#--plot growth matrices + mean growth (not an assignment)----
plst = rCompTCMs::compareResults.Pop.MeanGrowthPlusProbs(list(tcsam=lstT));
plst$male
plst$male + ggplot2::scale_y_continuous(limits=c(NA,185))

#--plot survey biomass, total catch mortality, and exploitation rates on same plot----
#--NOTE: making two plots because scales are different

##--extract NMFS EBS survey biomass data for males----
dfrSB = rTCSAM02::getMDFR.Data.FleetTimeSeries(lstT,fleet.type="survey",data.type="biomass",catch.type="index") |> 
          dplyr::filter(fleet %in% c("NMFS M"),x=="male",m=="all maturity",s=="all shell");
##--extract observer total catch data----
dfrCB = rTCSAM02::getMDFR.Data.FleetTimeSeries(lstT,fleet.type="fishery",data.type="biomass",catch.type="total") |> 
          dplyr::filter(fleet %in% c("TCF","SCF","RKF","GF All"),x=="all sex",m=="all maturity",s=="all shell",y>=1991) |> 
          dplyr::group_by(y) |> dplyr::summarize(val=wtsUtilities::Sum(val)) |> dplyr::ungroup() |> 
          dplyr::right_join(tibble::tibble(y=1991:2024),by="y");
##--extract retained catch data----
dfrRB = rTCSAM02::getMDFR.Data.FleetTimeSeries(lstT,fleet.type="fishery",data.type="biomass",catch.type="retained") |> 
          dplyr::filter(fleet %in% c("TCF"),x=="all sex",m=="all maturity",s=="all shell") |> 
          dplyr::right_join(tibble::tibble(y=1965:2024),by="y");
##--plot above on same scale----
ggplot(dfrRB,aes(x=y,y=val)) + 
  geom_line() + 
  geom_line(data=dfrCB,colour="blue") + 
  geom_line(data=dfrSB,colour="green") + 
  labs(y="biomass (1000's t)") + theme(axis.title.x=element_blank())
##--can't really calculate an informative exploitation rate from above

##--extract model-estimated population biomass by x,m,z----
dfrPBxmz = rTCSAM02::getMDFR.Pop.Biomass(lstT,cast="x+m+z") |> dplyr::select(y,x,m,z,pop=val) |> dplyr::mutate(y=as.numeric(y),z=as.numeric(z));
##--extract model-estimated total catch mortality by x,m,z----
dfrCBxmz = rTCSAM02::getMDFR.Fisheries.CatchBiomass(lstT,category="total mortality",cast="x+m+z") |> dplyr::mutate(y=as.numeric(y),z=as.numeric(z)) |> 
             dplyr::group_by(y,x,m,z) |> dplyr::summarize(tm=wtsUtilities::Sum(val)) |> dplyr::ungroup();
##--combine above
dfrPCBxmz = dfrPBxmz |> dplyr::full_join(dfrCBxmz,by=c("y","x","m","z")) |> 
              dplyr::mutate(pop=ifelse(is.na(pop),0,pop),
                            tm=ifelse(is.na(tm),0,tm));
##--create function to plot pop bio, total catch mortality, and exploitation rate
p<-function(dfr){
  p1 = ggplot(dfr,aes(x=y)) + 
        geom_line(aes(y=pop),colour="blue") + 
        geom_line(aes(y=tm),colour="green") + 
        scale_x_continuous(breaks=seq(1965,2025,5),limits=c(1965,NA)) + 
        labs(y="biomass (1000's t)") + theme(axis.title.x=element_blank());
  p2 = ggplot(dfr,aes(x=y)) + 
        geom_line(aes(y=expl)) + 
        scale_x_continuous(breaks=seq(1970,2025,5),limits=c(1965,NA)) + 
        labs(y="exploitation rate") + theme(axis.title.x=element_blank());
  pg = cowplot::plot_grid(p1,p2,ncol=1,rel_heights=c(2,1))
  pg
}
#--plot all Tanner crab----
dfrPCB = dfrPCBxmz |> 
           dplyr::group_by(y) |> dplyr::summarize(pop=sum(pop),tm=sum(tm)) |> dplyr::ungroup() |> 
           dplyr::mutate(expl=tm/pop);
p(dfrPCB);
#--plot preferred males----
dfrPCB = dfrPCBxmz |> dplyr::filter(x=="male",z>=125) |> 
           dplyr::group_by(y) |> dplyr::summarize(pop=sum(pop),tm=sum(tm)) |> dplyr::ungroup() |> 
           dplyr::mutate(expl=tm/pop);
p(dfrPCB);
#--plot mature males----
dfrPCB = dfrPCBxmz |> dplyr::filter(x=="male",m=="mature") |> 
           dplyr::group_by(y) |> dplyr::summarize(pop=sum(pop),tm=sum(tm)) |> dplyr::ungroup() |> 
           dplyr::mutate(expl=tm/pop);
p(dfrPCB);








