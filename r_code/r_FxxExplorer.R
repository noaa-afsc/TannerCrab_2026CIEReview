#--Code to do a little down-and-dirty exploring regarding reference points for Tanner crab
##--NOTE: this is NOT the way the assessment approaches it and is, at best, a very poor 
##        approximation to the assessment and at worst simply wrong.
require(ggplot2);
require(rlang);
require(rCompTCMs);

#--get TCSAM02 model results
lstT = wtsUtilities::getObj("../model_files/TCSAM02/Results.RData");


#--extract prM2M 
dfrPrM2M = rCompTCMs::extractMDFR.Pop.PrM2M(list(tcsam=lstT)) |> 
             dplyr::filter(x=="male");

#--extract population-level model-predicted new shell male abundance-at-size
##--calculate predicted fraction mature vs. all new shell, 
##--compare with assessment-estimated prM2M curve
dfrPopAbd = rCompTCMs::extractMDFR.Pop.Abundance(list(tcsam=lstT),cast="y+x+m+s+z") |> 
              dplyr::filter(x=="male",s=="new shell",as.numeric(y)>2012) |> dplyr::select(y,m,z,val) |> 
              tidyr::pivot_wider(names_from="m",values_from="val") |> 
              dplyr::mutate(val=mature/(mature+immature),
                            z=as.numeric(z),
                            y=as.character(y));
ggplot(dfrPopAbd,aes(x=z,y=val,colour=y)) + 
  geom_line() + 
  geom_line(data=dfrPrM2M,linetype=3) + 
  geom_hline(yintercept=0,linetype=3) + 
  labs(x="size (mm CW)",y="value",color="year")

#--extract population-level mature male abundance-at-size
##--look at abundance of mature males over time vs. average
dfrPopMMA = rCompTCMs::extractMDFR.Pop.Abundance(list(tcsam=lstT),cast="y+x+m+z") |> 
              dplyr::filter(x=="male",m=="mature",as.numeric(y)>2012) |> dplyr::select(y,z,val) |> 
              dplyr::mutate(z=as.numeric(z),
                            y=as.character(y));
dfrAvgPopMMA = dfrPopMMA |> dplyr::group_by(z) |> dplyr::summarize(MMA=mean(val)) |> dplyr::ungroup();
ggplot(dfrPopMMA,aes(x=z,y=val,colour=y)) + 
  geom_line() + 
  geom_line(aes(y=MMA),data=dfrAvgPopMMA,color="black",linewidth=1) + 
  geom_hline(yintercept=0,linetype=3) + 
  labs(x="size (mm CW)",y="MMA (millions)",color="year")


#--extract model-predicted population-level mature male biomass-at-size (includes old shell as well as new shell)
##--calculate cumulative sum across size by year, total by year, and 
##--forward and reverse cumulative fractions by year.
dfrPopMMB = rCompTCMs::extractMDFR.Pop.Biomass(list(tcsam=lstT),cast="y+x+m+s+z") |> 
              dplyr::filter(x=="male",m=="mature",as.numeric(y)>2012) |> dplyr::select(y,z,val) |> 
              dplyr::mutate(z=as.numeric(z),
                            y=as.character(y)) |> 
              dplyr::arrange(y,z) |> 
              dplyr::group_by(y,z) |> 
              dplyr::summarize(val=sum(val)) |>   #--sum over shell condition (although could eliminate "s" from `cast` above)
              dplyr::ungroup() |> 
              dplyr::group_by(y) |> 
              dplyr::mutate(cum=cumsum(val),
                            tot=sum(val)) |> 
              dplyr::ungroup() |> 
              dplyr::mutate(fpct=cum/tot,
                            rpct=1-fpct);
#--calculate average quantities
dfrAvgPopMMB = dfrPopMMB |> dplyr::group_by(z) |> 
                 dplyr::summarize(MMB=mean(val),
                                  cumMMB=mean(cum),
                                  fpct=mean(fpct),
                                  rpct=mean(rpct)) |> 
                 dplyr::ungroup();
#--plot average forward, reverse cumulative fractions (blue), prM2M (black), 
##--and choice of annual fractions by size
p = function(dfrPopMMB,val1,val2=NULL){
    val1 = rlang::enquo(val1);
    val2 = rlang::enquo(val2);
    p = ggplot(dfrPopMMB,aes(x=z,colour=y)) + geom_line(aes(y=!!val1));
    if (rlang::is_symbol(val2)) p = p + geom_line(aes(y=!!val2),data=dfrPopMMB); #--if val2 was NULL, it is not a symbol after enquo'ing
    p = p + 
          geom_line(aes(y=rpct),data=dfrAvgPopMMB,linetype=1,color="blue", linewidth=1) + 
          geom_line(aes(y=fpct),data=dfrAvgPopMMB,linetype=1,color="blue", linewidth=1) + 
          geom_line(aes(y=val),data=dfrPrM2M,    linetype=1,color="black",linewidth=1) + 
          geom_hline(yintercept=c(0,0.35,0.5,0.65),linetype=3) + 
          geom_vline(xintercept=125,linetype=3) + 
          labs(x="size (mm CW)",y="value",color="year") + 
          wtsPlots::getStdTheme();
    return(p);
}
#p(dfrPopMMB,cum);  #--absolute scale for annual cumulative MMB doesn't work with prM2M and fractions
p(dfrPopMMB,fpct);  #--plot annual forward cumulative fractions MMB
p(dfrPopMMB,rpct);  #--plot annual reverse cumulative fractions MMB

ggplot(dfrPopMMB,aes(x=z,y=val,color=y)) + geom_line() +
    geom_line(aes(y=MMB),data=dfrAvgPopMMB,color="blue",linewidth=1) + 
    geom_line(aes(y=cumMMB),data=dfrAvgPopMMB,color="blue",linewidth=1) + 
    geom_vline(xintercept=125,linetype=3) + 
    labs(x="size (mm CW)",y="MMB (1,000's t)",color="year") + 
    wtsPlots::getStdTheme();

#--extract annual directed fishery (TCF) selectivity and retention functions
dscM = 0.321;
##--calculate "vulnerability" curve = sel*(ret * dscM*(1-ret)) : fraction fishing mortality by size (i.e., unscaled by fully-selected F)
dfrSels = rCompTCMs::extractMDFR.Fisheries.SelFcns(list(tcsam=lstT),fleets="TCF") |> 
             dplyr::filter(x=="male",as.numeric(y)>2012) |> dplyr::select(y,z,val) |> dplyr::mutate(type="sel");
dfrRets = rCompTCMs::extractMDFR.Fisheries.RetFcns(list(tcsam=lstT),fleets="TCF") |> 
             dplyr::filter(x=="male",as.numeric(y)>2012) |> dplyr::select(y,z,val) |> dplyr::mutate(type="ret");
dfrVul  = dplyr::bind_rows(dfrSels,dfrRets) |> tidyr::pivot_wider(names_from="type",values_from="val") |> 
             dplyr::mutate(vul=sel*(ret+(1-ret)*dscM)) |> tidyr::pivot_longer(c("sel","ret","vul"),names_to="type",values_to="val");
##--compare annual ret, sel, vul with prM2M and annual reverse cumulative fraction MMB
ggplot(dfrVul,aes(x=z,y=val,colour=type)) + 
  geom_line() + 
  #geom_line(aes(y=fpct),data=dfrPopMMB,color="blue",linewidth=1) + 
  geom_line(aes(y=rpct),data=dfrPopMMB,color="blue",linewidth=1) + 
  geom_line(data=dfrPrM2M,color="black",linewidth=1) + 
  facet_wrap(~y) + 
  labs(x="size (mm CW)",y="value",color="type") + 
  wtsPlots::getStdTheme();
##--annual vulnerability curves and the average reverse cumulative MMB fraction-by-size
ggplot(dfrVul |> dplyr::filter(type=="vul"),aes(x=z,y=val,colour=y)) + 
  geom_line() + 
  #geom_line(aes(y=fpct),data=dfrPopMMB,color="blue",linewidth=1) + 
  geom_line(aes(y=rpct),data=dfrAvgPopMMB,color="blue",linewidth=1) + 
  #geom_line(data=dfrPrM2M,color="black",linewidth=1) + 
  geom_vline(xintercept=125,linetype=3) + 
  geom_hline(yintercept=c(0,0.35,0.5,0.65,1),linetype=3) + 
  labs(x="size (mm CW)",y="value",color="year") + 
  wtsPlots::getStdTheme();

#--combine vulnerability curve and MMB----
##--calculate vulnerability scaled to total MMB, relative MMB,
dfrVB = (dfrPopMMB |> dplyr::inner_join(dfrVul |> tidyr::pivot_wider(names_from="type",values_from="val"),by=c("y","z"))) |> 
             dplyr::mutate(vMMB=vul*val) |> 
             dplyr::select(y,z,MMB=val,vMMB,vul) |> 
             dplyr::arrange(y,z) |>
             dplyr::group_by(y) |> 
             dplyr::mutate(cumMMB=cumsum(MMB),totMMB=sum(MMB),
                           cumVMMB=cumsum(vMMB),totVMMB=sum(vMMB)) |>
             dplyr::ungroup() |> 
             dplyr::mutate(sclVul=vul*totMMB,
                           relMMB=MMB/totMMB,
                           relVMMB=vMMB/totMMB,
                           fcumMMB=cumMMB/totMMB,  rcumMMB =1-fcumMMB,
                           fcumVMMB=cumVMMB/totMMB,rcumVMMB=rcumMMB);
##--pivot to long format
dfrVBL = dfrVB |> 
            tidyr::pivot_longer(c("totMMB","MMB","vMMB","relMMB","relVMMB","cumMMB","cumVMMB","fcumMMB","fcumVMMB","rcumMMB","rcumVMMB","sclVul","vul"),names_to="type",values_to="val");
##--calculate averages by column
dfrAvgVB = dfrVB |> dplyr::group_by(z) |> dplyr::summarize(dplyr::across(MMB:rcumVMMB,mean));
##--plot annual vulnerable MMB by size
###--cumulative MMB by size: dotted blue line
###--cumulative vulnerable MMB by size: solid blue line
###--NOTE: C_inf = max(cumMMB) - max(cumVulnerableMMB)
ggplot(dfrVB,aes(x=z,y=vMMB,color=y)) + geom_line() + 
    geom_line(aes(y=cumVMMB),data=dfrAvgVB,color="blue",linewidth=1) + 
    geom_line(aes(y=cumMMB), data=dfrAvgVB,color="blue",linewidth=1,linetype=3) + 
    geom_hline(yintercept=0,linetype=3) + 
    geom_vline(xintercept=125,linetype=3) + 
    labs(x="size (mm CW)",y="MMB (1,000's t)",color="year") + 
    wtsPlots::getStdTheme();
totAvgMMB = dfrAvgVB$totMMB[1];
B_100 = max(dfrAvgVB$cumMMB);   #--"unfished" equilibrium biomass (F's relatively small, so ok approx.?)
C_inf = max(dfrAvgVB$cumVMMB);  #--theoretical catch at infinite F (i.e., all vulnerable MMB)
minXX = (B_100-C_inf)/B_100;  #--minimum ratio (=XX in F_XX)    #--minimum attainable XX for F_XX given B_100 and C_inf
F_XX = function(XX){
         F = log(C_inf/(C_inf-B_100*(1-XX))); # = log((C_inf/B_100)/((C_inf/B_100)-(1-XX)) = -log(1-(1-XX)/r) where r = C_inf/B_100
         return(F);
      }
F_35 = F_XX(0.35); #--NA!! can't reach it!
F_60 = F_XX(0.60); #--1.79

p = function(dfr,scl) {
  p = ggplot(dfr,aes(x=z,y=val,color=type)) + 
        geom_line() + 
        geom_hline(yintercept=scl*c(0,0.35,0.5,0.65,1),linetype=3) + 
        geom_vline(xintercept=125,linetype=3) + 
        scale_y_continuous(name="MMB (kt)",sec.axis=sec_axis(~.*(1/totAvgMMB),name="scaled value")) + 
        labs(x="size (mm CW)")
        wtsPlots::getStdTheme();
  p;
}
p(dfrAvgVBL |> dplyr::filter(type %in% c("cumMMB","cumVMMB","sclVul")),scl=totAvgMMB)

#--equilibrium quantities----
dfrM = rTCSAM02::getMDFR("ptrOFLResults/popDyInfoM/M_msz",lstT,TRUE);

##--get model-estimated population abundance quantities at unfished equilibrium----
dfrNatZ_F0 = rTCSAM02::getMDFR("ptrOFLResults/eqNatZF0_xmsz",lstT,TRUE) |> 
               dplyr::filter(x=="male",m=="mature") |> 
               dplyr::group_by(z) |> 
               dplyr::summarize(MMA=sum(val)) |> 
               dplyr::ungroup() |> 
               dplyr::mutate(z=as.numeric(z)) |> dplyr::arrange(z);
##--get model-estimated population abundance quantities at FM equilibrium----
dfrNatZ_FM = rTCSAM02::getMDFR("ptrOFLResults/eqNatZFM_xmsz",lstT,TRUE) |> 
               dplyr::filter(x=="male",m=="mature") |> 
               dplyr::group_by(z) |> 
               dplyr::summarize(MMA=sum(val)) |> 
               dplyr::ungroup() |> 
               dplyr::mutate(z=as.numeric(z)) |> dplyr::arrange(z);
ggplot(dfrNatZ_F0,aes(x=z,y=MMA)) + 
  geom_line() + 
  geom_line(data=dfrNatZ_FM,color="red") + 
  geom_line(data=dfrAvgPopMMA,color="blue");
##--get weights-at-size
dfrWatZ = rTCSAM02::getMDFR.Pop.WatZ(lstT) |> dplyr::filter(x=="male",m=="mature") |> 
            dplyr::mutate(z=as.numeric(z)) |> dplyr::arrange(z) |> 
            dplyr::select(z,wAtZ=val);
ggplot(dfrWatZ,aes(x=z,y=wAtZ)) + geom_line();
##--combine to calculate equilibrium unfished MMB-at-size
dfrBatZ_F0 = dfrNatZ_F0 |> dplyr::inner_join(dfrWatZ,by="z") |> 
               dplyr::mutate(MMB=wAtZ*MMA) |> 
               dplyr::mutate(cumMMB=cumsum(MMB));
##--combine to calculate equilibrium MMB-at-size at FM
dfrBatZ_FM = dfrNatZ_FM |> dplyr::inner_join(dfrWatZ,by="z") |> 
               dplyr::mutate(MMB=wAtZ*MMA) |> 
               dplyr::mutate(cumMMB=cumsum(MMB));
ggplot(mapping=aes(x=z)) + 
  geom_line(aes(y=MMB),data=dfrBatZ_F0,color="black") + 
  geom_line(aes(y=MMB),data=dfrBatZ_FM,color="red") + 
  geom_line(aes(y=MMB),data=dfrAvgPopMMB,color="blue") + 
  geom_line(aes(y=cumMMB),data=dfrBatZ_F0,color="black") + 
  geom_line(aes(y=cumMMB),data=dfrBatZ_FM,color="red") + 
  geom_line(aes(y=cumMMB),data=dfrAvgPopMMB,color="blue");
max(dfrBatZ_F0$cumMMB);
max(dfrBatZ_FM$cumMMB);
max(dfrBatZ_FM$cumMMB)/max(dfrBatZ_F0$cumMMB);
0.35*max(dfrBatZ_F0$cumMMB);

##--combine with average vulnerability curve
dfrEqSel = rTCSAM02::getMDFR("ptrOFLResults/catchInfoM/selF_fmsz",lstT,TRUE) |> 
             dplyr::filter(fleet=="TCF",m=="mature",s=="new shell") |> 
             dplyr::select(z,sel=val) |> dplyr::mutate(z=as.numeric(z)) |> dplyr::arrange(z);
dfrEqRet = rTCSAM02::getMDFR("ptrOFLResults/catchInfoM/retF_fmsz",lstT,TRUE) |> 
             dplyr::filter(fleet=="TCF",m=="mature",s=="new shell") |> 
             dplyr::select(z,ret=val) |> dplyr::mutate(z=as.numeric(z)) |> dplyr::arrange(z);
dfrEqSRV = dplyr::bind_cols(dfrEqSel,dfrEqRet |> dplyr::select(ret)) |> 
             dplyr::mutate(vul=sel*(ret+(1-ret)*dscM));

dfrBatZ_F0 = dplyr::bind_cols(dfrBatZ_F0,dfrEqSRV |> dplyr::select(!z)) |> 
               dplyr::mutate(vMMB=vul*MMB) |> 
               dplyr::mutate(cumVB=cumsum(vMMB));
B_100e = max(dfrBatZ_F0$cumMMB); ##--actual equilibrium unfished MMB  (150.6)
C_infe = max(dfrBatZ_F0$cumVB);  ##--actual max possible catch of MMB (81.209) (was 78.3)
minXX = (B_100e-C_infe)/B_100e;  ##--minimum possible XX              (0.46)   (was 0.48)
ggplot(dfrBatZ_F0,aes(x=z)) + 
        geom_line(aes(y=cumMMB),color="blue") + 
        geom_line(aes(y=cumVB), color="green") + 
        geom_line(aes(y=cumMMB),data=dfrBatZ_FM,color="red") + 
        geom_hline(yintercept=B_100e*c(0,0.35,0.5,0.65,1),linetype=3) + 
        geom_hline(yintercept=0,linetype=3) + 
        geom_vline(xintercept=125,linetype=3) + 
        scale_y_continuous(name="MMB (kt)",sec.axis=sec_axis(~.*(1/B_100e),name="scaled value")) + 
        labs(x="size (mm CW)") + 
        wtsPlots::getStdTheme();
B_FM = max(dfrBatZ_FM$cumMMB);
B_FM/B_100e;                        #-- 0.50
B_FM/lstT$rep$ptrOFLResults$B100;   #-- 0.61


lstMRPs = lstT$rep$ptrOFLResults;
lstMRPs$B100; #--123.482
lstMRPs$Fmsy; #--1.469
lstMRPs$Bmsy; #--43.219
lstMRPs$avgRec; #--601.213
lstMRPs$ofl_fx;

#----results with lstTp----
##--get model-estimated population abundance quantities at unfished equilibrium----
dfrNatZ_F0m = rTCSAM02::getMDFR("ptrOFLResults/eqNatZF0m_xmsz",lstTp,TRUE) |> 
               dplyr::filter(x=="male",m=="mature") |> 
               dplyr::group_by(z) |> 
               dplyr::summarize(MMA=sum(val)) |> 
               dplyr::ungroup() |> 
               dplyr::mutate(z=as.numeric(z)) |> dplyr::arrange(z);
##--get model-estimated population abundance quantities at FM equilibrium----
dfrNatZ_FMm = rTCSAM02::getMDFR("ptrOFLResults/eqNatZFMm_xmsz",lstTp,TRUE) |> 
               dplyr::filter(x=="male",m=="mature") |> 
               dplyr::group_by(z) |> 
               dplyr::summarize(MMA=sum(val)) |> 
               dplyr::ungroup() |> 
               dplyr::mutate(z=as.numeric(z)) |> dplyr::arrange(z);
ggplot(mapping=aes(x=z,y=MMA)) + 
  geom_line(data=dfrNatZ_F0,color="black") + 
  geom_line(data=dfrNatZ_F0,color="black",linetype=3) + 
  geom_point(data=dfrNatZ_F0,color="black") + 
  geom_line(data=dfrNatZ_FMm,color="red") + 
  geom_line(data=dfrNatZ_FM,color="red",linetype=3) + 
  geom_point(data=dfrNatZ_FM,color="red") + 
  geom_line(data=dfrAvgPopMMA,color="blue");
##--get weights-at-size
dfrWatZp = rTCSAM02::getMDFR.Pop.WatZ(lstTp) |> dplyr::filter(x=="male",m=="mature") |> 
            dplyr::mutate(z=as.numeric(z)) |> dplyr::arrange(z) |> 
            dplyr::select(z,wAtZ=val);
ggplot(dfrWatZp,aes(x=z,y=wAtZ)) + geom_line();
##--combine to calculate equilibrium unfished MMB-at-size
dfrBatZ_F0m = dfrNatZ_F0m |> dplyr::inner_join(dfrWatZp,by="z") |> 
               dplyr::mutate(MMB=wAtZ*MMA) |> 
               dplyr::mutate(cumMMB=cumsum(MMB));
##--combine to calculate equilibrium MMB-at-size at FM
dfrBatZ_FMm = dfrNatZ_FMm |> dplyr::inner_join(dfrWatZp,by="z") |> 
               dplyr::mutate(MMB=wAtZ*MMA) |> 
               dplyr::mutate(cumMMB=cumsum(MMB));
ggplot(mapping=aes(x=z)) + 
  geom_line(aes(y=MMB),data=dfrBatZ_F0m,color="black") + 
  geom_line(aes(y=MMB),data=dfrBatZ_FMm,color="red") + 
  geom_line(aes(y=MMB),data=dfrAvgPopMMB,color="blue") + 
  geom_line(aes(y=cumMMB),data=dfrBatZ_F0m,color="black") + 
  geom_line(aes(y=cumMMB),data=dfrBatZ_FMm,color="red") + 
  geom_line(aes(y=cumMMB),data=dfrAvgPopMMB,color="blue");
max(dfrBatZ_F0m$cumMMB);                         #--123.4982
max(dfrBatZ_FMm$cumMMB);                         #--43.22436
max(dfrBatZ_FMm$cumMMB)/max(dfrBatZ_F0m$cumMMB); #--0.35
0.35*max(dfrBatZ_F0m$cumMMB);                    #--43.22436

##--combine with average vulnerability curve
dfrEqSelm = rTCSAM02::getMDFR("ptrOFLResults/catchInfoM/selF_fmsz",lstTp,TRUE) |> 
             dplyr::filter(fleet=="TCF",m=="mature",s=="new shell") |> 
             dplyr::select(z,sel=val) |> dplyr::mutate(z=as.numeric(z)) |> dplyr::arrange(z);
dfrEqRetm = rTCSAM02::getMDFR("ptrOFLResults/catchInfoM/retF_fmsz",lstTp,TRUE) |> 
             dplyr::filter(fleet=="TCF",m=="mature",s=="new shell") |> 
             dplyr::select(z,ret=val) |> dplyr::mutate(z=as.numeric(z)) |> dplyr::arrange(z);
dfrEqSRVm = dplyr::bind_cols(dfrEqSelm,dfrEqRetm |> dplyr::select(ret)) |> 
             dplyr::mutate(vul=sel*(ret+(1-ret)*dscM));

dfrBatZ_F0m = dplyr::bind_cols(dfrBatZ_F0m,dfrEqSRVm |> dplyr::select(!z)) |> 
               dplyr::mutate(vMMB=vul*MMB) |> 
               dplyr::mutate(cumVB=cumsum(vMMB));

B_100m = max(dfrBatZ_F0m$cumMMB); ##--actual equilibrium unfished MMB  (123.4982)
C_infm = max(dfrBatZ_F0m$cumVB);  ##--actual max possible catch of MMB (66.5739)
minXXm = (B_100m-C_infm)/B_100m;  ##--minimum possible XX              (0.46)  <--WRONG!
ggplot(dfrBatZ_F0m,aes(x=z)) + 
        geom_line(aes(y=cumMMB),color="blue") + 
        geom_line(aes(y=cumVB), color="green") + 
        geom_line(aes(y=cumMMB),data=dfrBatZ_FM,color="red") + 
        geom_hline(yintercept=B_100e*c(0,0.35,0.5,0.65,1),linetype=3) + 
        geom_hline(yintercept=0,linetype=3) + 
        geom_vline(xintercept=125,linetype=3) + 
        scale_y_continuous(name="MMB (kt)",sec.axis=sec_axis(~.*(1/B_100e),name="scaled value")) + 
        labs(x="size (mm CW)") + 
        wtsPlots::getStdTheme();
B_FM = max(dfrBatZ_FM$cumMMB);
B_FM/B_100e;                        #-- 0.50
B_FM/lstT$rep$ptrOFLResults$B100;   #-- 0.61







