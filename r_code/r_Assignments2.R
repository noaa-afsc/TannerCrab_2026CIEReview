#--CIE Review assignments: Day 2
##--"Active Project" should be 2026-05_TannerCrab
require(ggplot2);

#--load project info----
fn = file.path(rstudioapi::getActiveProject(),paste0("rda_project_info.",wtsUtilities::getOperatingSystem(),".RData"));
s = wtsUtilities::getObj(fn);

#--load 2025 assessment model----
lstT = wtsUtilities::getObj(s$results$Asmt2025)

#--plot survey biomass, total catch mortality, and exploitation rates on same plot----
#--NOTE: making two plots because scales are different

##--extract NMFS EBS survey biomass data for males----
dfrSB = rTCSAM02::getMDFR.Data.FleetTimeSeries(lstT,fleet.type="survey",data.type="biomass",catch.type="index") |> 
          dplyr::filter(fleet %in% c("NMFS M"),x=="male",m=="all maturity",s=="all shell") |> 
          dplyr::select(y,val,lci,uci);
##--extract total catch data----
dfrCB = rTCSAM02::getMDFR.Data.FleetTimeSeries(lstT,fleet.type="fishery",data.type="biomass",catch.type="total") |> 
          dplyr::filter(fleet %in% c("TCF","SCF","RKF","GF All"),x=="all sex",m=="all maturity",s=="all shell",y>=1991) |> 
          dplyr::group_by(y) |> dplyr::summarize(total=wtsUtilities::Sum(val)) |> dplyr::ungroup() |> 
          dplyr::right_join(tibble::tibble(y=1991:2024),by="y");
##--extract retained catch data----
dfrRB = rTCSAM02::getMDFR.Data.FleetTimeSeries(lstT,fleet.type="fishery",data.type="biomass",catch.type="retained") |> 
          dplyr::filter(fleet %in% c("TCF","SCF","RKF","GF All"),x=="all sex",m=="all maturity",s=="all shell") |> 
          dplyr::group_by(y) |> dplyr::summarize(retained=wtsUtilities::Sum(val)) |> dplyr::ungroup() |> 
          dplyr::right_join(tibble::tibble(y=1965:2024),by="y");
ggplot(mapping=aes(x=y)) + 
  geom_line(aes(y=retained),data=dfrRB,colour="black") + 
  geom_line(aes(y=total),   data=dfrCB,colour="blue") + 
  geom_line(aes(y=val),     data=dfrSB,colour="green") + 
  geom_ribbon(aes(y=val,ymin=lci,ymax=uci),data=dfrSB,fill="green",colour=NA,alpha=0.2) + 
  labs(y="biomass (1000's t)") + theme(axis.title.x=element_blank());
ymin=1990;
ggplot(mapping=aes(x=y)) + 
  geom_line(aes(y=retained),data=dfrRB,colour="black") + 
  geom_line(aes(y=total),   data=dfrCB,colour="blue") + 
  geom_line(aes(y=val),     data=dfrSB,colour="green") + 
  geom_ribbon(aes(y=val,ymin=lci,ymax=uci),data=dfrSB,fill="green",colour=NA,alpha=0.2) + 
  scale_x_continuous(limits=c(ymin,NA)) + 
  scale_y_continuous(limits=c(0,max((dfrSB$uci[dfrSB$y>=ymin])))) + 
  labs(y="biomass (1000's t)") + theme(axis.title.x=element_blank())
dfrCBp = dfrRB |> dplyr:: inner_join(dfrCB,by="y") |> dplyr::filter(!is.na(retained),!is.na(total));
medR = median(dfrCBp$retained);
rngR = range(dfrCBp$retained);
sdvR = sd(dfrCBp$retained);
medC = median(dfrCBp$total);
rngC = range(dfrCBp$total);
sdvC = sd(dfrCBp$total);
medS = median(dfrSB$val[dfrSB$y>=ymin]);
sclS = ((dfrSB$val)-medS)/medS;
tacR1 = medR + ((dfrSB$val[dfrSB$y>=ymin])-medS)/medS *(rngR[2]-rngR[1]);
tacR2 = medR + ((dfrSB$val[dfrSB$y>=ymin])-medS)/medS *sdvR;
tacC1 = medC + ((dfrSB$val[dfrSB$y>=ymin])-medS)/medS *(rngC[2]-rngC[1]);
tacC2 = medC + ((dfrSB$val[dfrSB$y>=ymin])-medS)/medS *sdvC;
dfrMgt = tibble::tibble(y=1991:2025,tacR1=tacR1,tacR2=tacR2,tacC1=tacC1,tacC2=tacC2) |> 
           tidyr::pivot_longer(2:5,names_to="type",values_to="value");
p = ggplot(mapping=aes(x=y)) + 
      geom_point(aes(y=value,colour=type),data=dfrMgt,size=2,alpha=1) + 
      geom_line(aes(y=retained),data=dfrRB,colour="black") + 
      geom_line(aes(y=total),   data=dfrCB,colour="blue") + 
      geom_line(aes(y=val),     data=dfrSB,colour="green") + 
      geom_ribbon(aes(y=val,ymin=lci,ymax=uci),data=dfrSB,fill="green",colour=NA,alpha=0.2) + 
      geom_hline(yintercept=medS) + 
      scale_x_continuous(limits=c(ymin,NA)) + 
      scale_y_continuous(limits=c(0,max((dfrSB$uci[dfrSB$y>=ymin])))) + 
      labs(y="biomass (1000's t)") + theme(axis.title.x=element_blank())
p;
p + scale_y_log10()



##--extract model-estimated recruitment and mature biomass----
scl = 10;
dfrR_SSB = rTCSAM02::getMDFR.SdRep.RecAndSSB(lstT) |> 
             dplyr::select(variable,y,x,est,lci,uci) |> 
             dplyr::mutate(est=ifelse(variable=="lnRec",exp(est)/scl,est),
                           lci=ifelse(variable=="lnRec",exp(lci)/scl,lci),
                           uci=ifelse(variable=="lnRec",exp(uci)/scl,uci),
                           variable=ifelse(variable=="lnRec","R",paste(x,"MB")));
ggplot(dfrR_SSB |> dplyr::filter(y>=1982),aes(x=y,y=est,ymin=lci,ymax=uci,colour=variable,fill=variable)) + 
  geom_ribbon(colour=NA,alpha=0.2) + geom_line() + 
  geom_hline(yintercept=0,linetype=3) + 
  labs(x="year",y=paste0(scl," millions or 1,000's t"));

dfrR_SSBp1 = dfrR_SSB |> dplyr::select(!c(lci,uci,x)) |> tidyr::pivot_wider(id_cols=y,names_from="variable",values_from=est);
lag=4;
dfrR = dfrR_SSBp1 |> dplyr::select(y,R) |> dplyr::mutate(y=y-lag);
dfrR_SSBp = dfrR_SSBp1 |> dplyr::select(!R) |> dplyr::inner_join(dfrR,by="y") |> 
              dplyr::mutate(RMMB=R/`male MB`,
                            RMFB=R/`female MB`,
                            MMBR=1/RMMB,             #--wrong!!
                            MFBR=1/RMFB) |>          #--wrong!!
              dplyr::filter(y>=1982);
p1 = ggplot(dfrR_SSBp,aes(x=y,y=RMMB)) + 
      geom_hline(yintercept=0,linetype=3) + 
      geom_point() + labs(y="R/MMB",x="year");
p2 = ggplot(dfrR_SSBp,aes(x=y,y=RMFB)) + 
      geom_hline(yintercept=0,linetype=3) + 
      geom_point() + labs(y="R/MFB",x="year");
p = cowplot::plot_grid(p1,p2,ncol=1)  
p  

p1 = ggplot(dfrR_SSBp,aes(x=`male MB`,y=RMMB,colour=y)) + 
      geom_hline(yintercept=0,linetype=3) + 
      geom_point() + labs(y="R/MMB");
p2 = ggplot(dfrR_SSBp,aes(x=`female MB`,y=RMFB,colour=y)) + 
      geom_hline(yintercept=0,linetype=3) + 
      geom_point() + labs(y="R/MFB");
p = cowplot::plot_grid(p1,p2,ncol=1)  
p  







