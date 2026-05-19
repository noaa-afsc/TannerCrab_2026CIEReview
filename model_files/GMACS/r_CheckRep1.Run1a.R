#--Compare gmacs model run with TCSAM02 model run
require(ggplot2)
require(rlang)
require(wtsUtilities); #--see https://github.com/wStockhausen/wtsUtilities
require(wtsGMACS);     #--see https://github.com/wStockhausen/wtsGMACS

#--set up folder paths----
##--first: set working directory to top-level gmacs folder
dirThs = getwd();

#--source functions----
source(file.path(dirThs,"r_Functions.R"));

#--read in TCSAM02 model results----
##--define path to TCSAM02 results (change to your path)
path_to_TCSAM02_Results.RData = "../TCSAM02/Results.RData";
lstT = wtsUtilities::getObj(path_to_TCSAM02_Results.RData);

#--read GMACS model results----
##--this assumes you've run gmacs with gmacs.rep1 in sub-folder run1
if (TRUE){
  #--onlny need to run this section once
  fnRep1 = file.path(dirThs,"run1","gmacs.rep1");
  iln=1;
  repG =  wtsGMACS::readGmacsRep1(fnRep1,verbose=FALSE);
  wtsUtilities::saveObj(repG,file.path(dirThs,"run1","rda_gmacs_rep.RData"));
} else {
  repG = wtsUtilities::getObj(file.path(dirThs,"run1","rda_gmacs_rep.RData"));
}
##--if you didn't run the model, a link to "rda_gmacs_rep.RData" is provided on the Review webpage
### repG = wtsUtilities::getObj(file.path(dirThs,"rda_gmacs_rep.RData"));
lstGs = list();
class(lstGs) = "gmacs_reslst";
lstGs[["repsLst"]] = list();
lstGs[["repsLst"]][["gmacs"]] = repG;
single = TRUE;
sym = names(lstGs$repsLst)[1]; #--name of first GMACS model
gmacs=rlang::data_sym(sym);    #--use "!!gmacs" in dplyr code to get value of "sym"

#--model configuration----
zBs = repG$size_midpoints;
nZBs = length(zBs);

#--estimated parameters----
dfrEstPars = repG$`Estimated parameters` |> 
               dplyr::mutate(
                 dplyr::across(!par_type,as.numeric)
               );
View(dfrEstPars);
##--check parameters at bounds----
dfrPsAtBs = dfrEstPars |> 
              dplyr::filter((!is.na(status)) & (status!=0));
View(dfrPsAtBs)

#--population processes----
##--allometry: OK!----
dfrG = lstGs$repsLst$gmacs$IndivWeightByXM |> 
         dplyr::mutate(case="gmacs",.before=1) |> 
         tidyr::pivot_longer(tidyselect::all_of(as.character(seq(27,182,5))),
                             names_to="z",values_to="val") |> 
         dplyr::select(case,y=Year,x=Sex,m=Maturity,z,val) |> 
         dplyr::mutate(dplyr::across(c(2,5,6),as.numeric)) |> 
         dplyr::filter(y==2000) |> dplyr::select(!y);
dfrT = rTCSAM02::getMDFR.Pop.WatZ(lstT) |> 
        dplyr::select(case,x,m,z,val) |> 
        dplyr::mutate(z=as.numeric(z));
dfr = dplyr::bind_rows(dfrG,dfrT);
p = function(dfr){
  ggplot(dfr,aes(x=z,y=val,colour=case,shape=case)) + 
      geom_line() + geom_point() + 
      facet_wrap(x~m) + 
      labs(x="size (mm CW)",y="weight (kg)") + 
  wtsPlots::getStdTheme();
}
p(dfr);
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       case="tcsam-gmacs",val=tcsam-gmacs,
                       pdif=2*val/(gmacs+tcsam));
p(dfrp);

##--natural mortality: OK!----
dfr = wtsGMACS::extractNaturalMortality(lstGs,list(tcsam=lstT));
dfr = dfr |> dplyr::mutate(y=ifelse(y<1980,"normal",ifelse(y>1984,"normal","heightened"))) |> 
         dplyr::distinct(case,y,x,m,s,val) |> 
         dplyr::arrange(y,x,m,s,case);
p = function(dfr,includeRef=TRUE){
  p = ggplot(dfr,mapping=aes(x=x,y=val,colour=case,fill=case)) + 
    geom_col(position="dodge") + 
    geom_hline(yintercept=0.0,linetype=3) + 
    #scale_y_continuous(limits=c(0,mx)) + 
    facet_wrap(~y,ncol=1) + 
    labs(x="Population class",y="M") + 
    wtsPlots::getStdTheme() + 
    theme(axis.title.x=element_blank(),
         legend.direction="horizontal",
         legend.position="inside",
         legend.position.inside=c(0.01,0.99),
         legend.justification.inside=c(0,1));
  if (includeRef) p = p + geom_hline(yintercept=0.23,linetype=3);
  p
}
p(dfr |> dplyr::filter(m=="immature"));
p(dfr |> dplyr::filter(m=="mature"));
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=gmacs-tcsam,
                       pdif=100*2*val/(gmacs+tcsam));
p(dfrp |> dplyr::filter(m=="immature"),includeRef=FALSE);
p(dfrp |> dplyr::filter(m=="mature"),includeRef=FALSE);

##--growth----
###--mean growth: OK----
dfr = extractMeanGrowth(lstGs,lstT);
dfrp = dfr |> dplyr::select(case,x,z,val) |>
         tidyr::pivot_wider(names_from="case",values_from="val") |>
         dplyr::mutate(case="tcsam-gmacs",
                       val=(!!gmacs)-tcsam);
p = function(dfr,diff=FALSE){
  p = ggplot(dfr,aes(x=z,y=val,colour=case,shape=case)) + 
        geom_line() + geom_point() + 
        labs(x="pre-molt size (mm CW)",y="mean post-molt size") + 
        facet_grid(x~.) + 
        wtsPlots::getStdTheme();
  if (diff){
    p = p + geom_hline(yintercept=0,linetype=3);
  } else  {
    p = p + geom_abline(slope=1,linetype=3);
  }
  p
}
p(dfrp,diff=TRUE);

###--growth matrices: OK!!----
dfr = extractGrowthMatrices(lstGs,lstT);
dfrp = dfr |> dplyr::select(case,x,z,zp,val) |>
         tidyr::pivot_wider(names_from="case",values_from="val") |>
         dplyr::mutate(case="tcsam-gmacs",
                       val=(!!gmacs)-tcsam,
                       pdif=100*2*val/((!!gmacs)+tcsam));
p = function(dfr,diff=FALSE){
  p = ggplot(dfr ,aes(x=zp,y=val,colour=case,shape=case)) + 
        geom_line() + geom_point() + 
        facet_wrap(~z) + 
        labs(x="post-molt size (mm CW)",y="growth probability\ngiven pre-molt size") + 
          wtsPlots::getStdTheme();
  if (diff) p = p + geom_hline(yintercept=0,linetype=3);
  p;
}
p(dfr |> dplyr::filter(x=="male"));
p(dfr |> dplyr::filter(x=="female"));
p(dfrp |> dplyr::filter(x=="male"),diff=TRUE);  #--OK!
p(dfrp |> dplyr::filter(x=="female"),diff=TRUE);#--OK!

##--probability of maturing: OK!!----
dfr = extractPrM2M(lstGs,lstT) |> 
        dplyr::filter(((y==2024)&(case=="gmacs")|(case=="tcsam"))) |>
        dplyr::select(case,x,z,val);
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=tcsam-gmacs,
                       pdif=100*2*val/(gmacs+tcsam));
p = function(dfr) {
  p = ggplot(dfr,mapping=aes(x=z,y=val,colour=case,shape=case)) + 
        geom_line() + geom_point() + 
    geom_hline(yintercept=0,linetype=3) + 
        labs(x="size (mm CW)",y="probability") + 
        wtsPlots::getStdTheme();
  p
}
p(dfr |> dplyr::filter(x=="male"));
p(dfr |> dplyr::filter(x=="female"));
p(dfrp |> dplyr::filter(x=="male"));
p(dfrp |> dplyr::filter(x=="female"));

#--fishery processes----
##--fishery retention curves: OK!----
dfrG = lstGs$repsLst$gmacs$selfcns |> dplyr::mutate(case="gmacs",.before=1) |> 
         dplyr::filter(fleet %in% c("TCF"),type=="retained",sex=="male") |> 
         dplyr::mutate(dplyr::across(5+1:length(seq(27,182,5)),as.numeric)) |> 
         dplyr::rowwise() |> 
         dplyr::mutate(keep=sum(dplyr::c_across(5+1:length(seq(27,182,5)))),.before=1) |> 
         dplyr::filter(keep>0) |> 
         tidyr::pivot_longer(tidyselect::all_of(as.character(seq(27,182,5))))  |> 
         dplyr::select(!c(keep,type)) |> 
         dplyr::select(case,y=year,fleet,x=sex,z=name,val=value) |> 
         dplyr::filter(((fleet=="TCF"))) |>
         dplyr::mutate(z=as.numeric(z),
                       val=as.numeric(val));
dfrT = rTCSAM02::getMDFR.Fisheries.RetFcns(lstT$rep);
dfr = dplyr::bind_rows(dfrG |> dplyr::select(case,y,fleet,x,z,val),
                       dfrT |> dplyr::select(case,y,fleet,x,z,val)) |> 
        dplyr::mutate(y=as.numeric(y)) |> 
        dplyr::distinct(case,y,fleet,x,z,val) |> 
        dplyr::arrange(case,y,fleet,x,z);
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=tcsam-gmacs,
                       pdif=100*2*val/(gmacs+tcsam));
p = function(dfr,diff=FALSE){
  p = ggplot(dfr,mapping=aes(x=z,y=val,colour=case,shape=case)) + 
        geom_line() + geom_point() + 
        facet_wrap(~y) + 
        geom_hline(yintercept=0,linetype=3) + 
        labs(x="size (mm CW)",y="Retention Curves",shape="year") + 
        wtsPlots::getStdTheme();
  if (!diff) p = p + geom_hline(yintercept=1,linetype=3);
  p;
}
p(dfr |> dplyr::filter(y %in% c(1965,1991,2005,2015,2020,2024)));
p(dfrp |> dplyr::filter(y %in% c(1965,1991,2005,2015,2020,2024)),diff=TRUE);
# p = function(dfr,y,ytitle) {
#   ys = rlang::ensym(y);
#   ggplot(dfr,aes(x=z,y=!!ys)) + 
#     geom_line() + geom_point() + 
#     facet_wrap(~y) + 
#     geom_hline(yintercept=0,linetype=3) + 
#     geom_vline(xintercept=125,linetype=3) + 
#     labs(x="size(mm CW)",y=ytitle) + 
#     wtsPlots::getStdTheme();
# }
# p(dfrp,pdif,"Percent difference in retention functions (gmacs-tcsam)");

##--fishery selectivity (capture) curves: OK!!----
dfrG = lstGs$repsLst$gmacs$selfcns |> dplyr::mutate(case="gmacs",.before=1) |> 
         dplyr::filter(!(fleet %in% c("NMFS","BSFRF")),type=="capture") |> 
         dplyr::mutate(fleet=ifelse(fleet=="GF_All","GF All",fleet)) |> 
         dplyr::mutate(dplyr::across(5+1:length(seq(27,182,5)),as.numeric)) |> 
         dplyr::rowwise() |> 
         dplyr::mutate(keep=sum(dplyr::c_across(5+1:length(seq(27,182,5)))),.before=1) |> 
         dplyr::filter(keep>0) |> 
         tidyr::pivot_longer(tidyselect::all_of(as.character(seq(27,182,5))))  |> 
         dplyr::select(!c(keep,type)) |> 
         dplyr::select(case,y=year,fleet,x=sex,z=name,val=value) |> 
         #dplyr::filter(((fleet=="TCF"))) |>
         dplyr::mutate(z=as.numeric(z),
                       val=as.numeric(val));
dfrT = rTCSAM02::getMDFR.Fisheries.SelFcns(lstT$rep);
dfr = dplyr::bind_rows(dfrG |> dplyr::select(case,y,fleet,x,z,val),
                       dfrT |> dplyr::select(case,y,fleet,x,z,val)) |> 
        dplyr::mutate(y=as.numeric(y)) |> 
        dplyr::distinct(case,y,fleet,x,z,val) |> 
        dplyr::arrange(case,y,fleet,x,z);
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=tcsam-gmacs,
                       pdif=100*2*val/(gmacs+tcsam));
p = function(dfr,diff=FALSE){
  p = ggplot(dfr,mapping=aes(x=z,y=val,colour=case,shape=case,group=paste(case,as.character(y)))) + 
      geom_line() + geom_point() + 
      geom_vline(xintercept=125,linetype=3) + 
      geom_hline(yintercept=0,linetype=3) + 
      labs(x="size (mm CW)",y="Fishery Selectivity",shape="year") + 
      wtsPlots::getStdTheme();
  if (!diff) p = p + geom_hline(yintercept=1,linetype=3);
  if (length(unique(dfr$x))>1){
    p = p + facet_grid(x~y,scales="free_y");
  } else {
    p = p + facet_wrap(~x+y);
  }
  p;
}
p(dfr  |> dplyr::filter(fleet=="TCF",y %in% c(1970,1980,1995)));
p(dfrp |> dplyr::filter(fleet=="TCF",y %in% c(1970,1980,1995)),diff=TRUE); # OK!!
p(dfr  |> dplyr::filter(fleet=="TCF",x=="male",y %in% c(2005:2009,2013:2015,2017:2018,2020:2024)));
p(dfrp |> dplyr::filter(fleet=="TCF",x=="male",y %in% c(2005:2009,2013:2015,2017:2018,2020:2024)),diff=TRUE); # OK!!
p(dfr  |> dplyr::filter(fleet=="SCF",y %in% c(1990,2000,2015)));               # OK!!
p(dfrp |> dplyr::filter(fleet=="SCF",y %in% c(1990,2000,2015)),diff=TRUE);
p(dfr  |> dplyr::filter(fleet=="RKF",y %in% c(1990,2000,2015)));
p(dfrp |> dplyr::filter(fleet=="RKF",y %in% c(1990,2000,2015)),diff=TRUE);     # OK!!
p(dfr  |> dplyr::filter(fleet=="GF All",y %in% c(1980,1995,2015)));
p(dfrp |> dplyr::filter(fleet=="GF All",y %in% c(1980,1995,2015)),diff=TRUE);  # OK!!

##--fishery Fs: OK!!----
dfr = extractFisheryFs(lstGs,lstT,season=2) |> 
        dplyr::filter(!(fleet %in% c("NMFS","BSFRF"))) |> 
        dplyr::mutate(fleet=ifelse(fleet=="GF_All","GF All",fleet)) |>
        dplyr::select(case,fleet,y,x,val);
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=tcsam-gmacs,
                       pdif=2*val/(gmacs+tcsam));
dfrp1 = dfrp |> 
         dplyr::filter(((fleet=="GF All")&(y>=1973))| #--parameter-estimated F's
                       ((fleet=="RKF")   &(y>=1990))|
                       ((fleet=="SCF")   &(y>=1990))|
                       ((fleet=="TCF")   &(y>=1965)));
dfrp2 = dfrp |> 
         dplyr::filter(((fleet=="RKF")   &(y<1990))|  #--effort-based extrapolated F's
                       ((fleet=="SCF")   &(y<1990)));
p = function(dfr,log=TRUE){
  p = ggplot(dfr,mapping=aes(x=y,y=val,colour=case,shape=case)) + 
      geom_line() + geom_point() + 
      facet_grid(fleet~x,scales="free_y") + 
      labs(x="year",y="Fully-selected Capture F's") + 
      wtsPlots::getStdTheme();
  if (log) {
    p = p + scale_y_log10();
  } else {
    p = p + geom_hline(yintercept=0,linetype=3);
  }
  p;
}
p(dfr  |> dplyr::filter(x=="male"));
p(dfrp1 |> dplyr::filter(x=="male"),log=FALSE);  #--OK!!
p(dfrp2 |> dplyr::filter(x=="male"),log=FALSE);  #--OK!!
p(dfr |> dplyr::filter(x=="female"));
p(dfrp1 |> dplyr::filter(x=="female"),log=FALSE);  #--OK!!
p(dfrp2 |> dplyr::filter(x=="female"),log=FALSE);  #--OK!!

#--survey processes----
##--NMFS survey selectivity curves: OK----
dfrG = lstGs$repsLst$gmacs$selfcns |> dplyr::mutate(case="gmacs",.before=1) |> 
         dplyr::filter(fleet %in% c("NMFS"),type=="capture") |> 
         dplyr::mutate(dplyr::across(5+1:length(seq(27,182,5)),as.numeric)) |> 
         dplyr::rowwise() |> 
         dplyr::mutate(keep=sum(dplyr::c_across(5+1:length(seq(27,182,5)))),.before=1) |> 
         dplyr::filter(keep>0) |> 
         tidyr::pivot_longer(tidyselect::all_of(as.character(seq(27,182,5))))  |> 
         dplyr::select(!c(keep,type)) |> 
         dplyr::select(case,y=year,fleet,x=sex,z=name,val=value) |> 
         dplyr::filter(((fleet=="NMFS")&(y>=1975))) |>
         dplyr::mutate(z=as.numeric(z),
                       val=as.numeric(val));
dfrT = rTCSAM02::getMDFR.Surveys.SelFcns(lstT$rep) |> 
         dplyr::mutate(fleet=ifelse(stringr::str_starts(fleet,"NMFS"),"NMFS","NA")) |>
         dplyr::filter(fleet!="NA");
dfr = dplyr::bind_rows(dfrG |> dplyr::select(case,y,fleet,x,z,val),
                       dfrT |> dplyr::select(case,y,fleet,x,z,val)) |> 
        dplyr::mutate(y=as.numeric(y)) |> 
        dplyr::distinct(case,y,fleet,x,z,val) |> 
        dplyr::arrange(case,y,fleet,x,z);
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=gmacs-tcsam,
                       pdif=2*val/(gmacs+tcsam));
p = function(dfr,diff=FALSE){
  p = ggplot(dfr,mapping=aes(x=z,y=val,colour=case,shape=as.character(y),group=paste(x,as.character(y)))) + 
      geom_line() + geom_point() + 
      facet_wrap(~x,ncol=1) + 
      geom_hline(yintercept=0,linetype=3) + 
      labs(x="size (mm CW)",y="NMFS Survey Selectivity",shape="year") + 
      wtsPlots::getStdTheme();
  if (!diff) p = p+geom_hline(yintercept=1,linetype=3);
  p;
}
p(dfr  |> dplyr::filter(y %in% c(1975,2000)));
p(dfrp |> dplyr::filter(y %in% c(1975,2000)),diff=TRUE);

##--BSFRF survey "selectivity" curves: OK----
###--curves are really "availability"
dfrG = lstGs$repsLst$gmacs$selfcns |> dplyr::mutate(case="gmacs",.before=1) |> 
         dplyr::filter(fleet %in% c("BSFRF"),type=="capture") |> 
         dplyr::mutate(dplyr::across(5+1:length(seq(27,182,5)),as.numeric)) |> 
         dplyr::rowwise() |> 
         dplyr::mutate(keep=sum(dplyr::c_across(5+1:length(seq(27,182,5)))),.before=1) |> 
         dplyr::filter(keep>1.0e-08) |> 
         tidyr::pivot_longer(tidyselect::all_of(as.character(seq(27,182,5))))  |> 
         dplyr::select(!c(keep,type)) |> 
         dplyr::select(case,y=year,fleet,x=sex,z=name,val=value) |> 
         #dplyr::filter(((fleet=="BSFRF")&(y>=1975))) |>
         dplyr::mutate(z=as.numeric(z),
                       val=as.numeric(val));
dfrT = rTCSAM02::getMDFR.Surveys.AvlFcns(lstT$rep) |> 
         dplyr::mutate(fleet=ifelse(stringr::str_starts(fleet,"SBS BSFRF"),"BSFRF","NA")) |>
         dplyr::filter(fleet!="NA");
dfr = dplyr::bind_rows(dfrG |> dplyr::select(case,y,fleet,x,z,val),
                       dfrT |> dplyr::select(case,y,fleet,x,z,val)) |> 
        dplyr::mutate(y=as.numeric(y)) |> 
        #dplyr::distinct(case,y,fleet,x,z,val) |> 
        dplyr::arrange(case,y,fleet,x,z);
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=gmacs-tcsam,
                       pdif=2*val/(gmacs+tcsam));
p = function(dfr,diff=FALSE){
  p = ggplot(dfr,mapping=aes(x=z,y=val,colour=case)) + 
      geom_line() + geom_point() + 
      facet_grid(y~x) + 
      geom_hline(yintercept=0,linetype=3) + 
      labs(x="size (mm CW)",y="BSFRF Survey Availability") + 
      wtsPlots::getStdTheme();
  if (!diff) p = p+geom_hline(yintercept=1,linetype=3);
  p;
}
p(dfr);
p(dfrp,diff=TRUE);

##--fully-selected survey catchabilities: OK----
###--for GMACS, values are duplicated by maturity for females, 
###--so need to keep only values for mature females
dfrG = lstGs$repsLst$gmacs$Index_fit_summary |> dplyr::mutate(case="gmacs") |> 
         dplyr::filter(maturity!="immature") |>
         dplyr::select(case,y=year,fleet,x=sex,val=q) |> 
         dplyr::mutate(val=as.numeric(val)); 
dfrT = rTCSAM02::getMDFR.Surveys.Catchability(lstT$rep) |> 
         dplyr::mutate(fleet=ifelse(stringr::str_starts(fleet,"NMFS"),"NMFS",
                                    ifelse(stringr::str_starts(fleet,"SBS BSFRF"),"BSFRF","NA"))) |>
         dplyr::filter((fleet!="NA")&(y!=2020));#--drop missing survey year also
dfr = dplyr::bind_rows(dfrG |> dplyr::select(case,y,fleet,x,val),
                       dfrT |> dplyr::select(case,y,fleet,x,val)) |> 
        dplyr::mutate(y=as.numeric(y));
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=tcsam-gmacs,
                       pdif=2*val/(gmacs+tcsam));
p = function(dfr,diff=FALSE){
  p = ggplot(dfr,mapping=aes(x=y,y=val,colour=case,shape=x)) + 
      geom_line() + geom_point() + 
      facet_wrap(~fleet,ncol=1) + 
      geom_hline(yintercept=0,linetype=3) + 
      labs(x="survey year",y="fully-selected Q",shape="sex") + 
      wtsPlots::getStdTheme();
  if (!diff) p=p+geom_hline(yintercept=1,linetype=3);
  p;
}
p(dfr  |> dplyr::filter(fleet=="NMFS"));
p(dfrp |> dplyr::filter(fleet=="NMFS"),diff=TRUE);
p(dfr  |> dplyr::filter(fleet=="BSFRF"));
p(dfrp |> dplyr::filter(fleet=="BSFRF"),diff=TRUE);

#--population quantities----
##--recruitment: OK----
dfr = extractRecruitment(lstGs,lstT);
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=tcsam-gmacs,
                       pdif=100*2*val/(gmacs+tcsam));
p = function(dfr){
  ggplot(dfr,mapping=aes(x=y,y=val,colour=case,shape=case)) + 
    geom_line() + geom_point() + 
    geom_hline(yintercept=0,linetype=3) +
    geom_vline(xintercept=c(1948,1953,1965,1971,1973,1980,1990,1996,2005,2020),linetype=3) +
    labs(y='Total Recruitment (millions)') + 
    wtsPlots::getStdTheme() + theme(axis.title.x=element_blank());
}
p(dfr);
p(dfrp);
  
##--size at recruitment: OK----
dfr = extractSizeAtRecruitment(lstGs,lstT);
dfrp = dfr |> dplyr::filter(x!="female",y!="1948-1973") |> dplyr::select(case,z,val) |> 
         tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       x="all",
                       val=tcsam-gmacs,
                       pdif=100*2*val/(gmacs+tcsam));
p = function(dfr){
  ggplot(dfr,mapping=aes(x=z,y=val,colour=x,linetype=case)) + 
      geom_line() + geom_point() + 
      geom_hline(yintercept=0,linetype=3) +
      labs(x="size(mm CW)",y="proportion") +
      wtsPlots::getStdTheme();
}
p(dfr); 
p(dfrp);

##--cohort progression----
###--aggregated abundance by x/m/s: OK!----
dfr = dplyr::bind_rows(extractCohortProgression(lstGs,lstT,cast="x+m+s",gmacsType=4) |> dplyr::filter(x=="female"),
                       extractCohortProgression(lstGs,lstT,cast="x+m+s",gmacsType=4) |> dplyr::filter(x=="male")) |> 
      dplyr::mutate(ms=paste0(m,"\n",s)) |> 
      dplyr::filter(!((m=="immature")&(s=="old shell")));
dfrp = dfr |> dplyr::select(case,x,ms,y,val) |>
       dplyr::arrange(x,ms,y,case,val) |>
       tidyr::pivot_wider(names_from="case",values_from="val") |>
       dplyr::mutate(case="tcsam-gmacs",
                     val=tcsam-(!!gmacs),
                     pdif=100*2*(val)/((!!gmacs)+tcsam),
                     ln_ratio=log((!!gmacs)/tcsam));
p = function(dfr){
  ggplot(dfr,mapping=aes(x=y,y=val,colour=case,shape=case)) + 
      geom_line() + geom_point() + 
      labs(x="year",y="abundance (millions)") + 
      geom_hline(yintercept=0,linetype=3) +
      facet_grid(ms~x,scales="free_y") + 
      wtsPlots::getStdTheme();
}
p(dfr);
p(dfrp);
###--abundance-at-size: OK!----
dfr = dplyr::bind_rows(extractCohortProgression(lstGs,lstT,cast="x+m+s+z") |> dplyr::filter(x=="female"),
                       extractCohortProgression(lstGs,lstT,cast="x+m+s+z") |> dplyr::filter(x=="male")) |> 
      dplyr::mutate(ms=paste0(m,"\n",s),
                    z=as.numeric(z)) |> 
      dplyr::filter(!((m=="immature")&(s=="old shell"))&(y<11));
dfrp = dfr |> dplyr::select(case,x,ms,y,z,val) |>
       dplyr::arrange(x,ms,y,case,z,val) |>
       tidyr::pivot_wider(names_from="case",values_from="val") |>
       dplyr::mutate(case="tcsam-gmacs",
                     val=tcsam-(!!gmacs),
                     pdif=100*2*(val)/((!!gmacs)+tcsam),
                     ln_ratio=log((!!gmacs)/tcsam));
p = function(dfr){
  ggplot(dfr,mapping=aes(x=z,y=val,colour=case,shape=ms)) + 
      geom_line() + geom_point() + 
      labs(x="year",y="abundance (millions)") + 
      geom_hline(yintercept=0,linetype=3) +
      facet_wrap(~y,scales="free_y") + 
      wtsPlots::getStdTheme() + theme(axis.title.x=element_blank());
}
p(dfr  |> dplyr::filter(x=="male",ms=="immature\nnew shell"));
p(dfrp |> dplyr::filter(x=="male",ms=="immature\nnew shell")); #--OK
p(dfr  |> dplyr::filter(x=="male",ms=="mature\nnew shell"));
p(dfrp |> dplyr::filter(x=="male",ms=="mature\nnew shell"));   #--OK
p(dfr  |> dplyr::filter(x=="male",ms=="mature\nold shell"));
p(dfrp |> dplyr::filter(x=="male",ms=="mature\nold shell"));   #--OK
p(dfr  |> dplyr::filter(x=="female",ms=="immature\nnew shell"));
p(dfrp |> dplyr::filter(x=="female",ms=="immature\nnew shell")); #--OK
p(dfr  |> dplyr::filter(x=="female",ms=="mature\nnew shell"));
p(dfrp |> dplyr::filter(x=="female",ms=="mature\nnew shell"));   #--OK
p(dfr  |> dplyr::filter(x=="female",ms=="mature\nold shell"));
p(dfrp |> dplyr::filter(x=="female",ms=="mature\nold shell"));   #--OK

##--aggregated pop abundance: OK!!----
dfr = dplyr::bind_rows(extractPopAbd(lstGs,lstT,cast="x+m+s") |> dplyr::filter(x=="female"),
                       extractPopAbd(lstGs,lstT,cast="x+m+s") |> dplyr::filter(x=="male")) |> 
      dplyr::mutate(ms=paste0(m,"\n",s)) |> 
      dplyr::filter(!((m=="immature")&(s=="old shell")));
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=tcsam-gmacs);
p = function(dfr) {
      xints = c(1948,1953,1965,1971,1980,1990,1996,2005,2020);
      rng = range(dfr$y);
      xints = xints[(rng[1]<=xints)&(xints<=rng[2])]
      ggplot(dfr,mapping=aes(x=y,y=val,colour=case,shape=case)) + 
        geom_line() + geom_point() + 
        geom_hline(yintercept=0,linetype=3) + 
        geom_vline(xintercept=xints,linetype=3) + 
        labs(x="year",y="abundance (millions)") + 
        facet_grid(ms~x,scales="free_y") + 
        wtsPlots::getStdTheme() + theme(axis.title.x=element_blank());
}
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,1948,1965))); #--OK 
p(dfrp |> dplyr::filter(x=="male",y<=1974));                     #--OK
p(dfrp |> dplyr::filter(x=="male",y>=1975));                     #--OK
p(dfrp |> dplyr::filter(x=="female",dplyr::between(y,1948,1965))); #--OK
p(dfrp |> dplyr::filter(x=="female",y<=1974));                     #--OK
p(dfrp |> dplyr::filter(x=="female",y>=1975));                     #--OK 

##--pop abundance size comps: OK!!----
dfr = dplyr::bind_rows(extractPopAbd(lstGs,lstT,cast="x+m+s+z") |> dplyr::filter(x=="female"),
                       extractPopAbd(lstGs,lstT,cast="x+m+s+z") |> dplyr::filter(x=="male")) |> 
      dplyr::mutate(ms=paste0(m,"\n",s),
                    z=as.numeric(z)) |> 
      dplyr::filter(!((m=="immature")&(s=="old shell")));
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="%(tcsam-gmacs)",
                       val=200*(tcsam-gmacs)/(tcsam+gmacs)) |> 
         dplyr::arrange(y,x,ms,as.numeric(z));
#View(dfrp |> dplyr::filter(x=="male")); 
p = function(dfr){
  p = ggplot(dfr,mapping=aes(x=z,y=val,colour=case,shape=ms)) + 
       geom_line() + geom_point() + 
       geom_hline(yintercept=0,linetype=3) + 
       facet_grid(y~ms,scales="free_y") + 
       labs(x="size (mm CW)",y="abundance") + 
       wtsPlots::getStdTheme();
  p;
}
p(dfr |> dplyr::filter(x=="male",dplyr::between(y,1948,1953))); #--OK
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,1948,1953)));#--OK
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,1954,1960)));#--OK
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,1961,1965)));#--ok
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,1966,1970)));#--OK
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,1971,1975)));#--OK
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,2006,2010)));#--OK
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,2021,2025)));#--OK

###--check MMOs (mature male ogives) from exported pop-abundance: OK!!----
dfrGMMO = dfr |> dplyr::filter(case=="gmacs",x=="male",s!="old shell") |> 
                 dplyr::select(!c(x,s,ms)) |> 
                 tidyr::pivot_wider(names_from="m",values_from="val");
dfrTMMO = dfr |> dplyr::filter(case=="tcsam",x=="male",s!="old shell") |> 
                 dplyr::select(!c(x,s,ms)) |> 
                 tidyr::pivot_wider(names_from="m",values_from="val");
dfr = dplyr::bind_rows(dfrGMMO,dfrTMMO) |> 
        dplyr::mutate(tot=immature+mature,
                      prM=mature/tot);
difs = function(dfr,d){
  ds = rlang::ensym(d);
  dfrp = dfr |> dplyr::select(case,y,z,!!ds) |> 
                tidyr::pivot_wider(names_from=case,values_from=!!ds) |> 
                dplyr::mutate(!!ds:=(tcsam-gmacs)/(tcsam+gmacs),
                              case="%(tcsam-gmacs)");
  return(dfrp);
}
p = function(dfr,y){
  ys=rlang::ensym(y);
  p = ggplot(dfr,mapping=aes(x=z,y=!!ys,colour=case,fill=case)) +
      geom_line() + geom_point() + 
      geom_hline(yintercept=0,linetype=3) + 
      geom_vline(xintercept=125,linetype=3) + 
      facet_wrap(~y,scales="free_y") +
      labs(x="size (mm CW)",y=paste("male maturity",rlang::as_string(ys))) +
      wtsPlots::getStdTheme();
  p;
}
p(dfr |> dplyr::filter(dplyr::between(y,2006,2010)),mature);
p(difs(dfr |> dplyr::filter(dplyr::between(y,2006,2010)),mature),mature)
p(dfr |> dplyr::filter(dplyr::between(y,2006,2010)),prM);
p(difs(dfr |> dplyr::filter(dplyr::between(y,2006,2010)),prM),prM)


##--aggregated pop biomass: OK!!----
dfr = dplyr::bind_rows(extractPopBio(lstGs,lstT,cast="x+m+s") |> dplyr::filter(x=="female"),
                       extractPopBio(lstGs,lstT,cast="x+m+s") |> dplyr::filter(x=="male")) |>
      dplyr::mutate(ms=paste0(m,"\n",s)) |> 
      dplyr::filter(!((m=="immature")&(s=="old shell")));
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="%(tcsam-gmacs)",
                       val=200*(tcsam-gmacs)/(tcsam+gmacs));
p = function(dfr){
      xints = c(1948,1953,1965,1971,1980,1990,1996,2005,2020);
      rng = range(dfr$y);
      xints = xints[(rng[1]<=xints)&(xints<=rng[2])]
      p = ggplot(dfr,mapping=aes(x=y,y=val,colour=case,shape=case)) + 
            geom_line() + geom_point() + 
            geom_hline(yintercept=0,linetype=3) + 
            geom_vline(xintercept=xints,linetype=3) + 
            labs(x="year",y="biomass (1,000's t)") + 
            facet_grid(ms~x,scales="free_y") + 
            wtsPlots::getStdTheme() + theme(axis.title.x=element_blank());
      p;
}
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,1948,1965))); #--OK 
p(dfrp |> dplyr::filter(x=="male",y<=1974));                     #--OK
p(dfrp |> dplyr::filter(x=="male",y>=1975));                     #--OK
p(dfrp |> dplyr::filter(x=="female",dplyr::between(y,1948,1965))); #--OK
p(dfrp |> dplyr::filter(x=="female",y<=1974));                     #--OK
p(dfrp |> dplyr::filter(x=="female",y>=1975));                     #--OK

#--Predicted Survey MMOs (model size bins):OK----
zBs = lstGs$repsLst$gmacs$size_midpoints;
nZBs = length(zBs);
dfrG = lstGs$repsLst$gmacs$PredictedSurveyMMOs |> 
        tidyr::pivot_longer(cols=2+1:nZBs,names_to="z",values_to="val") |> 
        dplyr::select(y=year,type,z,val) |> 
        dplyr::mutate(dplyr::across(c(y,z,val),as.numeric)) |>
        tidyr::pivot_wider(names_from="type",values_from="val") |> 
        dplyr::mutate(case="gmacs",.before=1); 
ys = (dfrG |> dplyr::group_by(case,y) |> dplyr::summarize(tot=sum(totN)) |> dplyr::filter(tot>0))[["y"]];
dfrT = rTCSAM02::getMDFR.Pop.MaturityOgives(lstT) |> 
         dplyr::select(case,y,z,prMat=val) |> 
         dplyr::mutate(dplyr::across(c(y,z),as.numeric));
dfr = dplyr::bind_rows(dfrG,dfrT) |> dplyr::filter(y %in% ys);
p = function(dfr,y){
  y = rlang::enquo(y);
  ggplot(dfr,mapping=aes(x=z,y=!!y,colour=case)) + 
     geom_line() + geom_point(size=1) + 
     geom_hline(yintercept=0,linetype=3) + 
     geom_vline(xintercept=125,linetype=3) + 
     facet_wrap(~y,ncol=3) + 
     labs(x="size (mm CW)",y="pr(mature)") + 
     wtsPlots::getStdTheme();
}
p(dfr,prMat);
diff=function(dfr,d){
  d = rlang::enquo(d);
  dfrp = dfr |> dplyr::select(c(case,y,z,!!d)) |>
           tidyr::pivot_wider(names_from="case",values_from=!!d) |> 
           dplyr::arrange(y,z) |> 
           dplyr::mutate(case="%(tcsam-gmacs)",
                         val=2*(tcsam-gmacs)/(gmacs+tcsam));
  dfrp;
}
p(dfrp<-diff(dfr |> dplyr::filter(y>=2006),prMat),val); #--OK

#--Predicted Survey MMOs (observed size bins): BAD----
zBs = lstGs$repsLst$gmacs$mmodZBs;
nZBs = length(zBs);
dfrG = lstGs$repsLst$gmacs$PredicteSurveyMMOsP |> 
        tidyr::pivot_longer(cols=2+1:nZBs,names_to="z",values_to="val") |> 
        dplyr::select(y=year,type,z,val) |> 
        dplyr::mutate(dplyr::across(c(y,z,val),as.numeric)) |>
        tidyr::pivot_wider(names_from="type",values_from="val") |> 
        dplyr::mutate(case="gmacs",.before=1); 
ys = (dfrG |> dplyr::group_by(case,y) |> dplyr::summarize(tot=sum(totN)) |> dplyr::filter(tot>0))[["y"]];
dfrT = rTCSAM02::getMDFR.Fits.MaturityOgiveData(lstT) |> 
         dplyr::select(case,y,type,z,val) |> 
         dplyr::mutate(dplyr::across(c(y,z,val),as.numeric)) |> 
         tidyr::pivot_wider(names_from="type",values_from="val") |> 
        dplyr::mutate(case="tcsam",.before=1); 
dfr = dplyr::bind_rows(dfrG |> dplyr::filter(y %in% ys) |> 
                         dplyr::select(case,y,z,N=totN,prMat),
                      dfrT |> dplyr::select(case,y,z,N=n,prMat=predicted));
p = function(dfr,y){
  y = rlang::enquo(y);
  ggplot(dfr,mapping=aes(x=z,y=!!y,colour=case)) + 
     geom_line() + geom_point(size=1) + 
     geom_hline(yintercept=0,linetype=3) + 
     geom_vline(xintercept=125,linetype=3) + 
     facet_wrap(~y,ncol=3) + 
     labs(x="size (mm CW)",y="pr(mature)") + 
     wtsPlots::getStdTheme();
}
p(dfr,prMat);
p(dfrp<-diff(dfr,prMat),val); #--OK

#--predicted data quantities----
##--predicted fishery capture abundance: OK!!----
dfr = extractFisheryCaptureAbundance(lstGs,lstT) |> 
        dplyr::filter(!(fleet %in% c("NMFS","BSFRF"))) |> 
        dplyr::mutate(fleet=ifelse(fleet=="GF_All","GF All",fleet));
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=200*(tcsam-!!gmacs)/(tcsam+!!gmacs));
p = function(dfr,diff=FALSE){
  if (diff){
    title="% difference, Total Fishery Capture Abundance";
  } else {
    title="Total Fishery Capture Abundance (millions)";
  }
  p = ggplot(dfr,mapping=aes(x=y,y=val,colour=case,shape=case)) + 
      geom_line() + geom_point() + 
      facet_grid(fleet~x,scale="free_y") + 
      labs(x="year",y=title) + 
      wtsPlots::getStdTheme();
  p;
}
p(dfr  |> dplyr::filter(x=="male"));
p(dfrp |> dplyr::filter(x=="male"),TRUE);
p(dfr  |> dplyr::filter(x=="female"));
p(dfrp |> dplyr::filter(x=="female"),TRUE)

##--predicted fishery capture biomass----
dfr = extractFisheryCaptureBiomass(lstGs,lstT) |> 
        dplyr::filter(!(fleet %in% c("NMFS","BSFRF"))) |> 
        dplyr::mutate(fleet=ifelse(fleet=="GF_All","GF All",fleet));
p = function(dfr){
  p = ggplot(dfr,mapping=aes(x=y,y=val,colour=case,shape=case)) + 
      geom_line() + geom_point() + scale_y_log10() + 
      facet_grid(fleet~.,scale="free_y") + 
      labs(x="year",y="Total Fishery Capture Biomass (1,000's t??)") + 
      wtsPlots::getStdTheme();
}
p(dfr |> dplyr::filter(x=="male"));
p(dfr |> dplyr::filter(x=="female"));

##--predicted total fishing mortality (biomass): OK!!----
dfr = extractTotalFishingMortality(lstGs,lstT);
dfrp =  dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=200*(tcsam-!!gmacs)/(tcsam+!!gmacs));
p = function(dfr,diff=FALSE){
  if (diff){
    title="% difference(Total Fishing Mortality)";
  } else {
    title="Total Fishing Mortality (1,000's t)";
  }
  p = ggplot(dfr,aes(x=y,y=val,colour=case,shape=case)) + 
        geom_line() + geom_point() + 
        labs(x="year",y=title) + 
        wtsPlots::getStdTheme();
  p;
}
p(dfr);
p(dfrp,TRUE);

##--predicted retained catch mortality (biomass)----
dfr = extractRetainedCatchMortality(lstGs,lstT);
dfrp =  dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=200*(tcsam-!!gmacs)/(tcsam+!!gmacs));
ggplot(dfr,aes(x=y,y=val,colour=case,shape=case)) + 
  geom_line() + geom_point() + 
  geom_vline(xintercept=c(1973,1980,1990,1996),linetype=3) + 
  labs(x="year",y="Retained Catch Mortality (1,0000's t)") + 
  wtsPlots::getStdTheme();
ggplot(dfr |> dplyr::filter(y>=1965),aes(x=y,y=val,colour=case,shape=case)) + 
  geom_line() + geom_point() + 
  scale_y_log10() + 
  labs(x="year",y="Retained Catch Mortality (1,0000's t)") + 
  wtsPlots::getStdTheme();

##--predicted fishery discard mortality----
dfr = extractDiscardCatchMortality(lstGs,lstT) |> 
        dplyr::filter(!(fleet %in% c("NMFS","BSFRF"))) |> 
        dplyr::mutate(fleet=ifelse(fleet=="GF_All","GF All",fleet));
dfrp =  dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=200*(tcsam-!!gmacs)/(tcsam+!!gmacs));
p = function(dfr){
      p = ggplot(dfr,mapping=aes(x=y,y=val,colour=case,shape=case)) + 
      geom_line() + geom_point() + 
      geom_hline(yintercept=0,linetype=3) + 
      facet_wrap(~fleet,ncol=1,scales="free_y") + 
      labs(x="year",y="Fishery Discard Mortality (1,000's t)") + 
      wtsPlots::getStdTheme();
}
p(dfr);
p + (dfr |> dplyr::filter(y<1991));
p + (dfr |> dplyr::filter(y>1990));
p + (dfr |> dplyr::filter(y>1990,fleet=="TCF"));
p = p + scale_y_log10();
p + dfr;
p + (dfr |> dplyr::filter(y<1991));
p + (dfr |> dplyr::filter(y>1990));

##--predicted survey quantities
p = function(dfr,diff=FALSE){
  if (diff){
    title = "% difference Predicted Survey Biomass";
  } else {
    title = "Predicted Survey Biomass (1,000's t)";
  }
  p = ggplot(dfr,aes(x=y,y=val,colour=case,shape=case)) + 
        geom_line() + geom_point() + 
        geom_hline(yintercept=0,linetype=3) + 
        facet_wrap(~x+m,ncol=1,scales="free_y") + 
        labs(x="year",y=title) + 
        wtsPlots::getStdTheme();
  p;
}
###--predicted survey biomass: NMFS: OK!!----
dfr = extractPredictedSurveyBiomass(lstGs,lstT,fleetG="NMFS",fleetT="NMFS");
dfrp =  dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=200*(tcsam-!!gmacs)/(tcsam+!!gmacs));
p(dfrp,TRUE);

##--predicted survey biomass: BSFRF: OK!!----
dfr = extractPredictedSurveyBiomass(lstGs,lstT,fleetG="BSFRF",fleetT="SBS BSFRF");
dfrp =  dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="tcsam-gmacs",
                       val=200*(tcsam-!!gmacs)/(tcsam+!!gmacs));
p(dfrp,TRUE);

##--predicted fishery size comps----
###--retained catch----

##--predicted survey abundance size comps----
dfr = dplyr::bind_rows(extractPopAbd(lstGs,lstT,cast="x+m+s+z") |> dplyr::filter(x=="female"),
                       extractPopAbd(lstGs,lstT,cast="x+m+s+z") |> dplyr::filter(x=="male")) |> 
      dplyr::mutate(ms=paste0(m,"\n",s),
                    z=as.numeric(z)) |> 
      dplyr::filter(!((m=="immature")&(s=="old shell")));
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(case="%(tcsam-gmacs)",
                       val=200*(tcsam-gmacs)/(tcsam+gmacs)) |> 
         dplyr::arrange(y,x,ms,as.numeric(z));
#View(dfrp |> dplyr::filter(x=="male")); 
p = function(dfr){
  p = ggplot(dfr,mapping=aes(x=z,y=val,colour=case,shape=ms)) + 
       geom_line() + geom_point() + 
       geom_hline(yintercept=0,linetype=3) + 
       facet_grid(y~ms,scales="free_y") + 
       labs(x="size (mm CW)",y="abundance") + 
       wtsPlots::getStdTheme();
  p;
}
p(dfr |> dplyr::filter(x=="male",dplyr::between(y,1948,1953))); #--OK
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,1948,1953)));#--OK
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,1954,1960)));#--OK
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,1961,1965)));#--ok
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,1966,1970)));#--OK
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,1971,1975)));#--OK
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,2006,2010)));#--OK
p(dfrp |> dplyr::filter(x=="male",dplyr::between(y,2021,2025)));#--OK

###--check MMOs from exported pop-abundance: OK!!----
dfrGMMO = dfr |> dplyr::filter(case=="gmacs",x=="male",s!="old shell") |> 
                 dplyr::select(!c(x,s,ms)) |> 
                 tidyr::pivot_wider(names_from="m",values_from="val");
dfrTMMO = dfr |> dplyr::filter(case=="tcsam",x=="male",s!="old shell") |> 
                 dplyr::select(!c(x,s,ms)) |> 
                 tidyr::pivot_wider(names_from="m",values_from="val");
dfr = dplyr::bind_rows(dfrGMMO,dfrTMMO) |> 
        dplyr::mutate(tot=immature+mature,
                      prM=mature/tot);
difs = function(dfr,d){
  ds = rlang::ensym(d);
  dfrp = dfr |> dplyr::select(case,y,z,!!ds) |> 
                tidyr::pivot_wider(names_from=case,values_from=!!ds) |> 
                dplyr::mutate(!!ds:=(tcsam-gmacs)/(tcsam+gmacs),
                              case="%(tcsam-gmacs)");
  return(dfrp);
}
p = function(dfr,y){
  ys=rlang::ensym(y);
  p = ggplot(dfr,mapping=aes(x=z,y=!!ys,colour=case,fill=case)) +
      geom_line() + geom_point() + 
      geom_hline(yintercept=0,linetype=3) + 
      geom_vline(xintercept=125,linetype=3) + 
      facet_wrap(~y,scales="free_y") +
      labs(x="size (mm CW)",y=paste("male maturity",rlang::as_string(ys))) +
      wtsPlots::getStdTheme();
  p;
}
p(dfr |> dplyr::filter(dplyr::between(y,2006,2010)),mature);
p(difs(dfr |> dplyr::filter(dplyr::between(y,2006,2010)),mature),mature)
p(dfr |> dplyr::filter(dplyr::between(y,2006,2010)),prM);
p(difs(dfr |> dplyr::filter(dplyr::between(y,2006,2010)),prM),prM)


#--fits to data----
##--load summary likelihood components----
dfrGs = lstGs$repsLst$gmacs$Likelihood_summary;

##--fishery catch biomass data likelihoods summary----
### these won't agree because the two approaches handle the 
### log-scale constant terms (0.5*log(2*pi)) differently
### gmacs includes them, TCSAM02 does not
dfrGCDs = dfrGs |> dplyr::filter((term=="Catch_data")&(type!="summary")&(units!="numbers")) |> 
            dplyr::mutate(type=paste(type,"catch"),
                          fleet=stringr::str_replace(fleet,"_"," "),
                          case="gmacs") |> 
            dplyr::select(case,fleet,type,units,nll,objfun=objfun_value) |> 
            dplyr::mutate(dplyr::across(c(nll,objfun),as.numeric));
dfrT = rTCSAM02::getMDFR.OFCs.FleetData(lstT,category="fisheries") |>
         dplyr::filter(data.type!="n.at.z") |>
         dplyr::filter(!((data.type=="abundance")&(fleet!="GF All"))) |>
         dplyr::filter(!((catch.type=="retained catch")&(x=="FEMALE"))) |>
         dplyr::select(case,fleet,type=catch.type,units=data.type,nll,objfun) |> 
         dplyr::mutate(case="tcsam");
dfr = dplyr::bind_rows(dfrGCDs,dfrT);
p = function(dfr){
  p = ggplot(dfr,mapping=aes(x=fleet,y=objfun,colour=case,fill=case)) + 
      geom_col(position="dodge") + 
      geom_hline(yintercept=0,linetype=3) + 
      facet_grid(type~.) + 
      labs(y="NLL (catch biomass)") + 
      wtsPlots::getStdTheme();
  p;
}
p(dfr);

##--fishery catch abundance data likelihoods summary----
### these won't agree because the two approaches handle the 
### log-scale constant terms (0.5*log(2*pi)) differently
### gmacs includes them, TCSAM02 does not
dfrGCDs = dfrGs |> dplyr::filter((term=="Catch_data")&(type!="summary")&(units=="numbers")) |> 
            dplyr::mutate(type=paste(type,"catch"),
                          fleet=stringr::str_replace(fleet,"_"," "),
                          case="gmacs") |> 
            dplyr::select(case,fleet,type,units,nll,objfun=objfun_value) |> 
            dplyr::mutate(dplyr::across(c(nll,objfun),as.numeric));
dfrT = rTCSAM02::getMDFR.OFCs.FleetData(lstT,category="fisheries") |>
         dplyr::filter(data.type!="n.at.z") |>
         dplyr::filter(((data.type=="abundance")&(fleet=="GF All"))) |>
         dplyr::filter(!((catch.type=="retained catch")&(x=="FEMALE"))) |>
         dplyr::select(case,fleet,type=catch.type,units=data.type,nll,objfun) |> 
         dplyr::mutate(case="tcsam");
dfr = dplyr::bind_rows(dfrGCDs,dfrT);
p = function(dfr){
  p = ggplot(dfr,mapping=aes(x=fleet,y=objfun,colour=case,fill=case)) + 
      geom_col(position="dodge") + 
      geom_hline(yintercept=0,linetype=3) + 
      facet_grid(type~.) + 
      labs(y="NLL (catch abundance)") + 
      wtsPlots::getStdTheme();
  p;
}
p(dfr);

##--aggregated time series----
diff = function(dfr,d){
  ds = rlang::ensym(d);
  dfrp = dfr |> dplyr::filter(observed!=0) |> dplyr::select(case,fleet,y,!!ds) |> 
                tidyr::pivot_wider(names_from=case,values_from=!!ds) |> 
                dplyr::mutate(!!ds:=tcsam-gmacs,case="tcsam-gmacs");
  return(dfrp);
}
###--fishery catch time series----
p0 = function(dfr,y,ytitle){
      ys = rlang::ensym(y);
      p = ggplot(dfr,mapping=aes(x=y,y=!!ys,colour=case,fill=case)) + 
          geom_line() + geom_point() + 
          geom_hline(yintercept=0,linetype=3) + 
          geom_vline(xintercept=c(1990,1996,2005,2013),linetype=3) + 
          facet_grid(fleet~.,scales="free_y") + 
          labs(x="year",y=paste(rlang::as_string(ys),ytitle)) + 
          wtsPlots::getStdTheme();
      p
}
####--retained catch biomass time series: OK!!----
##### calculate consistent terms (residual,zscore,nll) for comparisons
#### gmacs: rsd is log-scale residual; nll includes 0.5*log(2*pi) term
##### tcsam: calculate cv, residual, nll consistent with gmacs
dfrGCDs = lstGs$repsLst$gmacs$Catch_fit_summary |> 
            dplyr::filter(fleet=="TCF",type=="retained",units=="biomass") |> 
            dplyr::select(y=year,fleet,observed=obs,cv,predicted=prd,rsd,nll) |> 
            dplyr::mutate(across(c(y,observed,cv,predicted,rsd,nll),as.numeric),
                          zscore=rsd/sqrt(log(1+cv^2)),
                          case="gmacs") |> 
            dplyr::select(case,fleet,y,observed,cv,predicted,residual=rsd,zscore,nll);
dfrT = rTCSAM02::getMDFR.AllScores.Biomass(lstT,
                                           fleet.type="fishery",
                                           catch.type="retained") |> 
        tidyr::pivot_wider(names_from=type,values_from=val) |> 
        dplyr::filter(x!="female") |> 
        dplyr::mutate(cv=sqrt(exp(stdv^2)-1),
                      residual=log(observed)-log(predicted),
                      nll=-dnorm(residual,sd=stdv,log=TRUE),
                      case="tcsam") |> 
        dplyr::select(case,fleet,y,observed,cv,predicted,residual,zscore=`z-score`,nll);
dfr = dplyr::bind_rows(dfrGCDs,dfrT);
p = function(dfr,y){ys=rlang::ensym(y);p0(dfr,!!ys,"(retained catch biomass)")};
p(dfr,observed);
p(dfr,cv);
p(diff(dfr,observed),observed);
p(diff(dfr,cv),cv);
p(dfr,predicted)+geom_point(aes(y=observed),data=dfr |> dplyr::filter(case=="gmacs"),colour="black");
p(dfr,residual);
p(dfr,zscore);
p(dfr,nll);
p(diff(dfr,residual),residual);
p(diff(dfr,nll),nll);            #--OK!!

####--total fishery catch biomass time series data: OK!!----
##### calculate consistent terms (residual,zscore,nll) for comparisons
#### gmacs: rsd is log-scale residual; nll includes 0.5*log(2*pi) term
##### tcsam: calculate cv, residual, nll consistent with gmacs
dfrGCDs = lstGs$repsLst$gmacs$Catch_fit_summary |> 
            dplyr::filter(type=="total",units=="biomass") |> 
            dplyr::select(y=year,fleet,observed=obs,cv,predicted=prd,rsd,nll) |> 
            dplyr::mutate(across(c(y,observed,cv,predicted,rsd,nll),as.numeric),
                          fleet=stringr::str_replace(fleet,"_"," "),
                          stdv=sqrt(log(1+cv^2)),
                          zscore=rsd/stdv,
                          nllp=-dnorm(rsd,sd=stdv,log=TRUE),
                          case="gmacs") |> 
            dplyr::select(case,fleet,y,observed,cv,predicted,stdv,residual=rsd,zscore,nll,nllp);
dfrT = rTCSAM02::getMDFR.AllScores.Biomass(lstT,
                                           fleet.type="fishery",
                                           catch.type="total") |> 
        tidyr::pivot_wider(names_from=type,values_from=val) |> 
        dplyr::filter(x!="female") |> 
        dplyr::mutate(cv=sqrt(exp(stdv^2)-1),
                      residual=log(observed)-log(predicted),
                      nll=-dnorm(residual,sd=stdv,log=TRUE),
                      nllp=nll,
                      case="tcsam") |> 
        dplyr::select(case,fleet,y,observed,cv,predicted,stdv,residual,zscore=`z-score`,nll,nllp);
dfr = dplyr::bind_rows(dfrGCDs,dfrT);
p = function(dfr,y){ys=rlang::ensym(y);p0(dfr,!!ys,"(total catch biomass)")};
p(dfr,observed);
p(dfr,cv);
p(diff(dfr|> dplyr::filter(observed!=0),observed),observed)
p(diff(dfr|> dplyr::filter(observed!=0),cv),cv)
p(dfr,predicted)+geom_line(aes(y=observed),data=dfr |> dplyr::filter(case=="tcsam"),colour="black");
p(dfr |> dplyr::filter(observed!=0),residual);
p(dfr |> dplyr::filter(observed!=0),zscore);
p(dfr |> dplyr::filter(observed!=0),nll);
p(diff(dfr |> dplyr::filter(observed!=0),residual),residual);
p(diff(dfr |> dplyr::filter(observed!=0),nll),nll);
p(diff(dfr |> dplyr::filter(observed!=0),nllp),nllp);
p(diff(dfr |> dplyr::filter(observed!=0),stdv),stdv);
p(diff(dfr |> dplyr::filter(observed!=0),cv),cv);

p(dfr |> dplyr::filter(case=="gmacs"),predicted)+geom_point(aes(y=observed),data=dfr |> dplyr::filter(case=="tcsam"),colour="black");

p(diff(dfr,cv),cv);
p(diff(dfr,residual),residual);
p(diff(dfr |> dplyr::filter(fleet=="GF All"),nll),nll);

####--total catch abundance time series data: OK!!----
dfrGCDs = lstGs$repsLst$gmacs$Catch_fit_summary |> 
            dplyr::filter(type=="total",units=="numbers") |> 
            dplyr::select(y=year,fleet,observed=obs,cv,predicted=prd,rsd,nll) |> 
            dplyr::mutate(across(c(y,observed,cv,predicted,rsd,nll),as.numeric),
                          fleet=stringr::str_replace(fleet,"_"," "),
                          zscore=rsd/sqrt(log(1+cv^2)),
                          case="gmacs") |> 
            dplyr::select(case,fleet,y,observed,cv,predicted,residual=rsd,zscore,nll);
dfrT = rTCSAM02::getMDFR.AllScores.Abundance(lstT,
                                           fleet.type="fishery",
                                           catch.type="total") |> 
        tidyr::pivot_wider(names_from=type,values_from=val) |> 
        dplyr::filter(fleet=="GF All") |> 
        dplyr::mutate(cv=sqrt(exp(sdobs^2)-1),
                      residual=log(observed)-log(predicted),
                      nll=-dnorm(residual,sd=stdv,log=TRUE),
                      case="tcsam") |> 
        dplyr::select(case,fleet,y,observed,cv,predicted,residual,zscore=`z-score`,nll);
dfr = dplyr::bind_rows(dfrGCDs,dfrT);
p = function(dfr,y){ys=rlang::ensym(y);p0(dfr,!!ys,"(total catch abundance)")};
p(dfr,observed);
p(dfr,cv);
p(dfr,predicted)+geom_line(aes(y=observed),data=dfr |> dplyr::filter(case=="tcsam"),colour="black");
p(dfr,residual);
p(dfr,zscore);
p(dfr,nll);
p(diff(dfr,nll),nll);

###--
dfrT = rTCSAM02::getMDFR.Fits.EffortData(lstT)

####--index (survey) catch biomass summary: OK!!----
dfrG = lstGs$repsLst$gmacs$Likelihood_summary |> 
         dplyr::filter(term=="Index_data",type!="summary") |> 
         dplyr::select(fleet,x=sex,m=maturity,nll) |> 
         dplyr::mutate(nll=as.numeric(nll),
                       m=ifelse(x=="male","undetermined",m)) |>   #--"mature" male is really undetermined male
         dplyr::mutate(case="gmacs",.before=1);
dfrT = rTCSAM02::getMDFR.OFCs.FleetData(lstT,category="surveys") |>
         dplyr::filter(data.type!="n.at.z") |>
         dplyr::filter((data.type!="abundance")&
                       (!(fleet |> stringr::str_starts("SBS NMFS")))&
                       (rmse!=0)) |>
         dplyr::select(-c(s,category,y)) |> 
         dplyr::mutate(x=tolower(x),
                       m=ifelse(m=="ALL_MATURITY","undetermined",m) |> tolower(),
                       fleet=ifelse(stringr::str_starts(fleet,"NMFS"),"NMFS",fleet),
                       fleet=ifelse(stringr::str_starts(fleet,"SBS BSFRF"),"BSFRF",fleet),
                       case="tcsam") |> 
         dplyr::select(case,fleet,x,m,nll);
dfr = dplyr::bind_rows(dfrG,dfrT) |> dplyr::mutate(xm=paste(m,x))
p = function(dfr){
  p = ggplot(dfr,mapping=aes(x=case,y=nll,color=case,fill=case)) +
      geom_col() +
      geom_hline(yintercept=0,linetype=3) +
      facet_grid(fleet~xm,scales="free_y") +
      labs(y="NLL") +
      wtsPlots::getStdTheme();
  p;
}
p(dfr);
dfrp = dfr |> tidyr::pivot_wider(names_from=case,values_from=nll) |> dplyr::mutate(nll=tcsam-gmacs,case="tcsam-gmacs")
p(dfrp) + facet_grid(fleet~xm,scales="fixed");

###--index (survey) catch time series----
p1 = function(dfr,y,ytitle){
      ys = rlang::ensym(y);
      dfrp = dfr |> dplyr::mutate(xm=paste(x,m));
      p = ggplot(dfrp,mapping=aes(x=y,y=!!ys,colour=case,fill=case)) + 
          geom_line() + geom_point() + 
          geom_hline(yintercept=0,linetype=3) + 
          facet_grid(xm~fleet,scales="free_y") + 
          labs(x="year",y=paste(rlang::as_string(ys),ytitle)) + 
          wtsPlots::getStdTheme();
      p
}
difs = function(dfr,d){
  ds = rlang::ensym(d);
  dfrp = dfr |> dplyr::filter(observed!=0) |> dplyr::select(case,fleet,y,x,m,!!ds) |> 
                tidyr::pivot_wider(names_from=case,values_from=!!ds) |> 
                dplyr::mutate(!!ds:=tcsam-gmacs,case="tcsam-gmacs");
  return(dfrp);
}
####--index (survey) catch biomass time series data----
dfrGSDs = lstGs$repsLst$gmacs$Index_fit_summary |> dplyr::filter(units=="biomass") |> 
            dplyr::select(y=year,fleet,x=sex,m=maturity,
                          observed=obs,cv=actual_CV,predicted=prd,zscore=prsn_res) |> 
            dplyr::mutate(dplyr::across(c(y,observed,cv,predicted,zscore),as.numeric),
                          sdobs=sqrt(log(1+cv^2)),
                          nll=0.5*log(2*pi*sdobs^2) + (0.5*zscore^2),
                          residual=log(observed)-log(predicted)) |> 
            dplyr::mutate(case="gmacs",.before=1);
dfrT = rTCSAM02::getMDFR.AllScores.Biomass(lstT,
                                           fleet.type="survey",
                                           catch.type="index") |> 
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
dfr = dplyr::bind_rows(dfrGSDs,dfrT) |> dplyr::mutate(xm=paste(m,x));
# p = function(dfr,y){
#       ys = rlang::ensym(y);
#       dfrp = dfr |> dplyr::mutate(xm=paste(x,m));
#       p = ggplot(dfrp,mapping=aes(x=y,y=!!ys,colour=case,fill=case)) + 
#           geom_line() + geom_point() + 
#           geom_hline(yintercept=0,linetype=3) + 
#           facet_grid(xm~.,scales="free_y") + 
#           labs(x="year") + 
#           wtsPlots::getStdTheme();
#       p
#     }
p = function(dfr,y){ys=rlang::ensym(y);p1(dfr,!!ys,"(survey biomass)")};
p(dfr |> dplyr::filter(fleet=="NMFS"),observed);
p(dfr |> dplyr::filter(fleet=="NMFS"),cv);
p(dfr |> dplyr::filter(fleet=="NMFS"),predicted);
p(dfr |> dplyr::filter(fleet=="NMFS"),residual);
p(dfr |> dplyr::filter(fleet=="NMFS"),zscore);
p(dfr |> dplyr::filter(fleet=="NMFS"),nll);
p(difs(dfr |> dplyr::filter(fleet=="NMFS"),nll),nll)

p(dfr |> dplyr::filter(fleet=="BSFRF"),observed);
p(dfr |> dplyr::filter(fleet=="BSFRF"),cv);
p(dfr |> dplyr::filter(fleet=="BSFRF"),predicted);
p(dfr |> dplyr::filter(fleet=="BSFRF"),residual);
p(dfr |> dplyr::filter(fleet=="BSFRF"),zscore);
p(dfr |> dplyr::filter(fleet=="BSFRF"),nll);
p(difs(dfr |> dplyr::filter(fleet=="BSFRF"),nll),nll)

###--fishery size composition data----
dfrGZDs = dfrGs |> dplyr::filter(term=="Size_data") |>
            dplyr::filter(!fleet %in% c("NMFS","BSFRF"),type!="summary") |>
            dplyr::mutate(case="gmacs",s="undetermined",nll.type="unknown") |>
            dplyr::select(!c(term,emphasis)) |>
            dplyr::select(case,fleet,catch.type=type,nll.type,x=sex,m=maturity,s,nll,objfun=objfun_value) |>
            dplyr::mutate(dplyr::across(c(nll,objfun),as.numeric),
                          fleet=ifelse(fleet=="GF_All","GF All",fleet),
                                       catch.type=paste(catch.type,"catch")) |>
            dplyr::group_by(case,fleet,catch.type) |>  #--summarize over x,m,s
            dplyr::summarize(nll=sum(nll),
                             objfun=sum(objfun)) |>
            dplyr::ungroup();
dfrT = rTCSAM02::getMDFR.OFCs.FleetData(lstT,category="fisheries") |>
         dplyr::filter((data.type=="n.at.z")) |>
         dplyr::group_by(fleet,catch.type) |>              #--summarize over x,m,s
         dplyr::summarize(nll=sum(nll,na.rm=TRUE),
                          objfun=sum(objfun,ma.rm=TRUE)) |>
         dplyr::ungroup() |>
         dplyr::mutate(case="tcsam",.before=1);
dfr = dplyr::bind_rows(dfrGZDs,dfrT);
p = function(dfr){
  p = ggplot(dfr,mapping=aes(x=fleet,y=objfun,colour=case,fill=case)) +
      geom_col(position="dodge") +
      facet_grid(catch.type~.,scales="free_y") +
      labs(x="sex",y="survey size comps NLLs") +
      wtsPlots::getStdTheme();
  p;
}
p(dfr);

###--survey size composition data----
dfrGZDs = dfrGs |> dplyr::filter(term=="Size_data") |>
            dplyr::filter(fleet %in% c("NMFS","BSFRF")) |>
            dplyr::mutate(case="gmacs",s="undetermined",nll.type="unknown") |>
            dplyr::select(!c(term,type,emphasis)) |>
            dplyr::select(case,fleet,nll.type,x=sex,m=maturity,s,nll,objfun=objfun_value) |>
            dplyr::mutate(dplyr::across(c(nll,objfun),as.numeric));
dfrT = rTCSAM02::getMDFR.OFCs.FleetData(lstT,category="surveys") |>
         dplyr::filter((data.type=="n.at.z") &
                      (!(fleet |> stringr::str_starts("SBS NMFS")))) |>
         dplyr::group_by(fleet,fit.type,nll.type,x,m,s) |>
         dplyr::summarize(nll=sum(nll,na.rm=TRUE),
                          objfun=sum(objfun,ma.rm=TRUE)) |>
         dplyr::ungroup() |>
         dplyr::mutate(x=tolower(x),
                       m=ifelse(m=="ALL_MATURITY","undetermined",tolower(m)),
                       s=ifelse(s=="ALL_SHELL","undetermined",tolower(s)),
                       fleet=ifelse(fleet |> stringr::str_starts("NMFS"),"NMFS",fleet),
                       fleet=ifelse(fleet |> stringr::str_starts("SBS"),"BSFRF",fleet)) |>
         dplyr::mutate(case="tcsam",.before=1);
dfr = dplyr::bind_rows(dfrGZDs,dfrT |> dplyr::select(!fit.type)) |>
        dplyr::mutate(xm=paste(x,m));
dfrp = dfr |> dplyr::select(!c(nll.type,nll)) |> 
         tidyr::pivot_wider(names_from=case,values_from=objfun) |> 
         dplyr::mutate(diff=tcsam-gmacs,case="tcsam-gmacs");
p = function(dfr,y) {
      y = rlang::enquo(y);
      p = ggplot(dfr,mapping=aes(x=x,y=!!y,colour=case,fill=case)) +
            geom_col(position="dodge") + 
            geom_hline(yintercept=0,linetype=3) + 
            facet_grid(fleet~xm,scales="free_y") +
            labs(x="sex",y="survey size comps NLLs") +
            wtsPlots::getStdTheme();
      p;
}
p(dfr,objfun);
p(dfrp,diff);

###--growth data: OK!!----
dfrGGDs = dfrGs |> dplyr::filter(term=="Growth_data",type=="growth_data") |>
            dplyr::filter(type!="summary") |>
            dplyr::select(x=sex,nll,wgt=emphasis,objfun=objfun_value) |>
            dplyr::mutate(case="gmacs",.before=1) |>
            dplyr::mutate(dplyr::across(c(nll,wgt,objfun),as.numeric));
dfrT = rTCSAM02::getMDFR.OFCs.GrowthData(lstT) |>
         dplyr::select(x,nll,wgt,objfun) |>
         dplyr::mutate(case="tcsam",.before=1);
dfr = dplyr::bind_rows(dfrGGDs,dfrT);
p = function(dfr,y){
  y = rlang::enquo(y);
  p = ggplot(dfr,mapping=aes(x=x,y=!!y,colour=case,fill=case)) +
        geom_col(position="dodge") + 
        geom_hline(yintercept=0,linetype=3) + 
        labs(x="sex",y="growth data NLLs") +
        wtsPlots::getStdTheme();
  p
}

p(dfr,objfun);
dfrp = dfr |> dplyr::select(!c(nll,wgt)) |> 
         tidyr::pivot_wider(names_from=case,values_from=objfun) |> 
         dplyr::mutate(diff=tcsam-gmacs,case="tcsam-gmacs");
p(dfrp,diff);

###--mature male ogives data----
dfrGMMOs = dfrGs |> dplyr::filter(term=="MMOD_data",type=="mmod_data") |> 
            dplyr::mutate(dplyr::across(c(nll,emphasis,objfun_value),as.numeric)) |> 
            dplyr::filter(type!="summary",nll>0) |>
            dplyr::select(y=year,nll,wgt=emphasis,objfun=objfun_value) |>
            dplyr::mutate(case="gmacs",.before=1) |>
            dplyr::mutate(dplyr::across(c(y,nll,wgt,objfun),as.numeric),
                          objfun=ifelse(objfun==0,NA,objfun));
dfrT = rTCSAM02::getMDFR.OFCs.MaturityOgiveData(lstT) |>
         dplyr::select(y,nll,wgt,objfun) |>
         dplyr::mutate(case="tcsam",.before=1);
minY = min(dfrT$y);
dfr = dplyr::bind_rows(dfrGMMOs,dfrT) |> dplyr::filter(y>=minY);
p = function(dfr,y){
  y=rlang::enquo(y);
  p = ggplot(dfr,mapping=aes(x=y,y=!!y,colour=case,fill=case)) +
      geom_line() + geom_point() + 
      labs(x="year",y="male maturity ogive data NLLs") +
      wtsPlots::getStdTheme();
  p
}
p(dfr,objfun);
dfrp = dfr |> dplyr::select(!c(nll,wgt)) |> 
         tidyr::pivot_wider(names_from=case,values_from=objfun) |> 
         dplyr::mutate(diff=tcsam-gmacs,case="tcsam-gmacs");
p(dfrp,diff);

###--mature male ogives data details: OK!!----
dfrG = lstGs$repsLst$gmacs$likeSuumaryMMOsP |> 
          dplyr::select(y=year,z=zb,n=ss,obs,prd,res,zscr,nll) |>
            dplyr::mutate(dplyr::across(c(y,z,n,obs,prd,res,zscr,nll),as.numeric)) |> 
            dplyr::mutate(case="gmacs",.before=1);
dfrT = rTCSAM02::getMDFR.Fits.MaturityOgiveData(lstT) |> 
         dplyr::select(y,type,z,val)|>
         dplyr::mutate(case="tcsam",.before=1) |> 
         tidyr::pivot_wider(names_from="type",values_from="val") |> 
         dplyr::select(case,y,z,n,obs=observed,prd=predicted,zscr=zscores,nll=nlls) |> 
         dplyr::mutate(res=obs-prd);
dfr = dplyr::bind_rows(dfrG,dfrT);
difs = function(dfr,d){
  ds = rlang::ensym(d);
  dfrp = dfr |> dplyr::select(case,y,z,!!ds) |> 
                tidyr::pivot_wider(names_from=case,values_from=!!ds) |> 
                dplyr::mutate(!!ds:=tcsam-gmacs,case="tcsam-gmacs");
  return(dfrp);
}
p = function(dfr,y){
  ys=rlang::ensym(y);
  p = ggplot(dfr,mapping=aes(x=z,y=!!ys,colour=case,fill=case)) +
      geom_line() + geom_point() + 
      geom_hline(yintercept=0,linetype=3) + 
      geom_vline(xintercept=125,linetype=3) + 
      facet_wrap(~y,scales="free_y") +
      labs(x="size (mm CW)",y=paste("male maturity ogive data",rlang::as_string(ys))) +
      wtsPlots::getStdTheme();
  p;
}
p(dfr,obs);
p(difs(dfr,obs),obs)
p(dfr,n);
p(difs(dfr,n),n)
p(dfr,prd);
p(difs(dfr,prd),prd)
p(dfr,res);
p(difs(dfr,res),res)
p(dfr,nll);
p(difs(dfr,nll),nll)

#--Reference Points----
dfrT = rTCSAM02::getMDFR.ManagementQuantities(lstT) |> 
         dplyr::mutate(case="tcsam") |> 
         dplyr::select(case,type,val)  |> 
         tidyr::pivot_wider(names_from="type",values_from="val") |> 
         dplyr::mutate(Bmsy=0.35*B100) |> 
         tidyr::pivot_longer(cols=c(-1),names_to="type",values_to="val");
nflt = lstGs$repsLst$gmacs$Number_of_fleets;
dfrG = lstGs$repsLst$gmacs$Derived_quatities |> 
         dplyr::select(type=param,val=est) |> 
         dplyr::mutate(dplyr::across(c(val),as.numeric)) |> 
         dplyr::filter(dplyr::row_number() %in% c(1:7,(7+nflt))) |> 
         tidyr::pivot_wider(names_from="type",values_from="val") |> 
         dplyr::mutate(case=sym,
                       avgRec=male_spr_rbar+female_spr_rbar,
                       OFL=`OFL(tot)`,
                       Fofl=`Fofl(1)`,
                       Fmsy=`Fmsy(1)`,
                       curB=NA,
                       prjB=`Bcurr/BMSY`*BMSY,
                       MSY=NA,
                       B100=male_spr_rbar*`SSSB/R(F=0)`/1000) |>  #--think SSB/R(F=0) is for males; units are grams, need KT so x 1000 if rec in millions
         dplyr::select(case,OFL,Fofl,prjB,curB,Fmsy,Bmsy=BMSY,MSY,B100,avgRec) |> 
         tidyr::pivot_longer(cols=c(-1),names_to="type",values_to="val");
dfr = dplyr::bind_rows(dfrG,dfrT) |> dplyr::filter(type!="B100");
dfrp = dfr |> tidyr::pivot_wider(names_from="case",values_from="val") |> 
         dplyr::mutate(`%(tcsam-gmacs)`=200*(tcsam-gmacs)/(tcsam+gmacs)) |> 
         dplyr::select(!c(gmacs,tcsam)) |> 
         tidyr::pivot_longer(cols=c(-1),names_to="case",values_to="val");
p = function(dfr){
  p = ggplot(dfr,aes(x=case,y=val,color=case,fill=case)) + 
        geom_col() + 
        geom_hline(yintercept=0,linetype=3) + 
        scale_y_continuous(expand=expansion(0.15,0)) + 
        facet_grid(type~.,scales="free_y") + 
        labs(y="value") + 
        wtsPlots::getStdTheme();
  p;
}
legPos = theme(legend.position="bottom")
p(dfr) + legPos;
p(dfrp)+ legPos;
  
  
  
  