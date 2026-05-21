#--make plots
require(ggplot2);
require(rlang);
require(rCompTCMs);

#--get TCSAM02 model results
lstT = wtsUtilities::getObj("../model_files/TCSAM02/Results.RData");

#--extract prM2M 
dfrPrM2M = rCompTCMs::extractMDFR.Pop.PrM2M(list(tcsam=lstT)) |> 
             dplyr::filter(x=="male");

#--extract population-level new shell male abundance-at-size
##--look at fraction mature vs. all new shell
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
##--look at fraction mature vs. all new shell
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


#--extract population-level mature male biomass-at-size
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
dfrAvgPopMMB = dfrPopMMB |> dplyr::group_by(z) |> 
                 dplyr::summarize(MMB=mean(val),
                                  cumMMB=mean(cum),
                                  fpct=mean(fpct),
                                  rpct=mean(rpct)) |> 
                 dplyr::ungroup();
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
#p(dfrPopMMB,cum);
p(dfrPopMMB,rpct);
ggplot(dfrPopMMB,aes(x=z,y=val,color=y)) + geom_line() +
    geom_line(aes(y=MMB),data=dfrAvgPopMMB,color="blue",linewidth=1) + 
    geom_line(aes(y=cumMMB),data=dfrAvgPopMMB,color="blue",linewidth=1) + 
    geom_vline(xintercept=125,linetype=3) + 
    labs(x="size (mm CW)",y="MMB (1,000's t)",color="year") + 
    wtsPlots::getStdTheme();

#--extract TCF selectivity and retention functions
dfrSels = rCompTCMs::extractMDFR.Fisheries.SelFcns(list(tcsam=lstT),fleets="TCF") |> 
             dplyr::filter(x=="male",as.numeric(y)>2012) |> dplyr::select(y,z,val) |> dplyr::mutate(type="sel");
dfrRets = rCompTCMs::extractMDFR.Fisheries.RetFcns(list(tcsam=lstT),fleets="TCF") |> 
             dplyr::filter(x=="male",as.numeric(y)>2012) |> dplyr::select(y,z,val) |> dplyr::mutate(type="ret");
dfrVul  = dplyr::bind_rows(dfrSels,dfrRets) |> tidyr::pivot_wider(names_from="type",values_from="val") |> 
             dplyr::mutate(vul=sel*(ret+(1-ret)*0.321)) |> tidyr::pivot_longer(c("sel","ret","vul"),names_to="type",values_to="val");
ggplot(dfrVul,aes(x=z,y=val,colour=type)) + 
  geom_line() + 
  #geom_line(aes(y=fpct),data=dfrPopMMB,color="blue",linewidth=1) + 
  geom_line(aes(y=rpct),data=dfrPopMMB,color="blue",linewidth=1) + 
  geom_line(data=dfrPrM2M,color="black",linewidth=1) + 
  facet_wrap(~y) + 
  labs(x="size (mm CW)",y="value",color="type") + 
  wtsPlots::getStdTheme();
ggplot(dfrVul |> dplyr::filter(type=="vul"),aes(x=z,y=val,colour=y)) + 
  geom_line() + 
  #geom_line(aes(y=fpct),data=dfrPopMMB,color="blue",linewidth=1) + 
  geom_line(aes(y=rpct),data=dfrAvgPopMMB,color="blue",linewidth=1) + 
  #geom_line(data=dfrPrM2M,color="black",linewidth=1) + 
  geom_vline(xintercept=125,linetype=3) + 
  geom_hline(yintercept=c(0,0.35,0.5,0.65,1),linetype=3) + 
  labs(x="size (mm CW)",y="value",color="year") + 
  wtsPlots::getStdTheme();

#--combine vulnerability curve and MMB
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
                           fcumVMMB=cumVMMB/totMMB,rcumVMMB=rcumMMB)# |> 
             #tidyr::pivot_longer(c("totMMB","MMB","vMMB","relMMB","relVMMB","cumMMB","cumVMMB","fcumMMB","fcumVMMB","rcumMMB","rcumVMMB","sclVul","vul"),names_to="type",values_to="val");
dfrAvgVB = dfrVB |> dplyr::group_by(z) |> dplyr::summarize(dplyr::across(MMB:rcumVMMB,mean));
ggplot(dfrVB,aes(x=z,y=vMMB,color=y)) + geom_line() + 
    geom_line(aes(y=cumVMMB),data=dfrAvgVB,color="blue",linewidth=1) + 
    geom_line(aes(y=cumMMB), data=dfrAvgVB,color="blue",linewidth=1,linetype=3) + 
    geom_vline(xintercept=125,linetype=3) + 
    labs(x="size (mm CW)",y="MMB (1,000's t)",color="year") + 
    wtsPlots::getStdTheme();
totAvgMMB = (dfrAvgVB |> dplyr::filter(type=="totMMB"))$val[1];
B_100 = max((dfrAvgVB |> dplyr::filter(type=="cumMMB"))$val);
C_inf = max((dfrAvgVB |> dplyr::filter(type=="cumVMMB"))$val);
minX = (B_100-C_inf)/B_100;
F_x = function(x){
         F = log(C_inf/(C_inf-B_100*(1-x)));
         return(F);
      }
F_35 = F_x(0.35)

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
p(dfrAvgVB |> dplyr::filter(type %in% c("cumMMB","cumVMMB","sclVul")),scl=totAvgMMB)


#--get numbers-at-size for unfished population
dfrNatZ_F0 = rTCSAM02::getMDFR("ptrOFLResults/eqNatZF0_xmsz",lstT,TRUE) |> 
               dplyr::filter(x=="male",m=="mature") |> 
               dplyr::group_by(z) |> 
               dplyr::summarize(MMA=sum(val)) |> 
               dplyr::ungroup() |> 
               dplyr::mutate(z=as.numeric(z)) |> dplyr::arrange(z);
ggplot(dfrNatZ_F0,aes(x=z,y=MMA)) + geom_line() + geom_line(data=dfrAvgPopMMA,color="blue");

dfrWatZ = rTCSAM02::getMDFR.Pop.WatZ(lstT) |> dplyr::filter(x=="male",m=="mature") |> 
            dplyr::mutate(z=as.numeric(z)) |> dplyr::arrange(z) |> 
            dplyr::select(z,wAtZ=val);
ggplot(dfrWatZ,aes(x=z,y=wAtZ)) + geom_line();
dfrBatZ_F0 = dfrNatZ_F0 |> dplyr::inner_join(dfrWatZ,by="z") |> 
               dplyr::mutate(MMB=wAtZ*MMA);
ggplot(dfrBatZ_F0,aes(x=z,y=MMB)) + geom_line() +;

