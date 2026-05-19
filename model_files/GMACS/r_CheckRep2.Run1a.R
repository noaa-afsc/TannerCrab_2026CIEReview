#--Compare gmacs model run with TCSAM02 model run
require(ggplot2)
require(rlang);
require(wtsUtilities); #--see https://github.com/wStockhausen/wtsUtilities
require(wtsGMACS);     #--see https://github.com/wStockhausen/wtsGMACS

#--set up folder paths----
##--first: set working directory to top-level gmacs folder
dirThs = getwd();

#--source functions----
source(file.path(dirThs,"r_Functions2.R"));

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


#--plot fits to index time series data----
dfr = getFitsToIndexTimeSeries(lstGs,lstT) 
p1 = plotFitsToIndexTSs(dfr,f="NMFS");      print(p1);
p2 = p1 + scale_y_log10();                  print(p2);
dfrp = dfr |> dplyr::select(case,y,fleet,x,m,xm,predicted) |> 
              tidyr::pivot_wider(names_from="case",values_from="predicted") |>
              dplyr::mutate(case="%(tcsam-gmacs)",
                            diff=tcsam-gmacs,
                            pdif=200*(diff)/(tcsam+gmacs));
ggplot(dfrp |> dplyr::filter(fleet=="NMFS"),aes(x=y,y=pdif,color=case)) + 
  geom_line() + 
  geom_hline(yintercept=0,linetype=3) + 
  facet_grid(xm~.) + 
  labs(y="% difference(predicted)") + 
  wtsPlots::getStdTheme() + wtsPlots::noXT();

#--plot fits to ZCs----
dfr = getFitsToZCs(lstGs,lstT);
plotFitsToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x=="male",  m=="undetermined",s=="undetermined"))
plotFitsToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x=="female",m=="mature",      s=="undetermined"))
plotFitsToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x=="female",m=="immature",    s=="undetermined"))
plotFitsToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="retained",x=="male",  m=="undetermined",s=="undetermined"))
plotFitsToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
plotFitsToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
plotFitsToZCs(dfr |> dplyr::filter(fleet=="SCF", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
plotFitsToZCs(dfr |> dplyr::filter(fleet=="SCF", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
plotFitsToZCs(dfr |> dplyr::filter(fleet=="RKF", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
plotFitsToZCs(dfr |> dplyr::filter(fleet=="RKF", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
plotFitsToZCs(dfr |> dplyr::filter(fleet=="GF All", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
plotFitsToZCs(dfr |> dplyr::filter(fleet=="GF All", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))

#--plot residuals to ZC fits----
plotResidualsToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x==  "male",m=="undetermined",s=="undetermined"))
plotResidualsToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x=="female",m=="mature",      s=="undetermined"))
plotResidualsToZCs(dfr |> dplyr::filter(fleet=="NMFS",comp_type=="total",   x=="female",m=="immature",    s=="undetermined"))
plotResidualsToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="retained",x=="male",  m=="undetermined",s=="undetermined"))
plotResidualsToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
plotResidualsToZCs(dfr |> dplyr::filter(fleet=="TCF", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
plotResidualsToZCs(dfr |> dplyr::filter(fleet=="SCF", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
plotResidualsToZCs(dfr |> dplyr::filter(fleet=="SCF", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
plotResidualsToZCs(dfr |> dplyr::filter(fleet=="RKF", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
plotResidualsToZCs(dfr |> dplyr::filter(fleet=="RKF", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))
plotResidualsToZCs(dfr |> dplyr::filter(fleet=="GF_All", comp_type=="total"   ,x=="male",  m=="undetermined",s=="undetermined"))
plotResidualsToZCs(dfr |> dplyr::filter(fleet=="GF_All", comp_type=="total",   x=="female",m=="undetermined",s=="undetermined"))

