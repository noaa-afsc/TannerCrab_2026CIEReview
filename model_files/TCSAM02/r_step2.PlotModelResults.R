
#-load model
require(rCompTCMs);
require(rTCSAM02);
resDir = "./best_results";
models<-list();
case<-c("2025");
models[[case]] = wtsUtilities::getObj(file.path(resDir,"Results.RData"));
case<-c("2024");
models[[case]] = wtsUtilities::getObj(file.path("../TCSAM02_2024_22_03d5/best_results","Results.RData"));
compare<-names(models);

#--extract tables and create plots 
dfr.PVs   <-rCompTCMs::extractMDFR.Results.ParameterValues(models[compare]);
readr::write_csv(dfr.PVs,file.path(resDir,"ParamValues.csv"))
dfr.PsAtBs<-rCompTCMs::extractMDFR.Results.ParametersAtBounds(models[compare]);
if (!is.null(dfr.PsAtBs)&&inherits(dfr.PsAtBs,"data.frame")){
  readr::write_csv(dfr.PsAtBs,file.path(resDir,"ParamsAtBounds.csv"));
  cat(paste0(dfr.PsAtBs$name," ",dfr.PsAtBs$test," [",dfr.PsAtBs$label,"]\n"));
} else {
  message("No parameters at bounds!!");
}

dfr.OFLs<-rTCSAM02::getMDFR.OFLResults(models[compare],verbose=TRUE);
if (!is.null(dfr.OFLs)){
  dfr.OFLs$avgRec<-dfr.OFLs$avgRecM + dfr.OFLs$avgRecF;
  dfr.OFLs<-wtsUtilities::deleteCols(dfr.OFLs,cols=c("avgRecM","avgRecF"),debug=FALSE);
  dfr.OFLs<-dfr.OFLs[,c("case","objFun","maxGrad","avgRec","B100","Bmsy","curB","Fmsy","MSY","Fofl","OFL","prjB")];
  readr::write_csv(dfr.OFLs,file.path(resDir,"oflResults.csv"))
}

rCompTCMs::compareResults.Surveys.SelFcns(models[compare],singlePlot=TRUE,
                                          dodge=0,years=c(1981,2018));
rCompTCMs::compareResults.Surveys.CaptureProbs(models[compare],singlePlot=TRUE,
                                               dodge=0,years=c(1981,2018));
rCompTCMs::compareResults.Fisheries.SelFcns(models[compare],fleets="RKF",singlePlot=TRUE,
                                            dodge=0,years=c(1990,2000,2015))
rCompTCMs::compareResults.Fisheries.SelFcns(models[compare],fleets="SCF",singlePlot=TRUE,
                                            dodge=0,years=c(1990,2000,2015))
rCompTCMs::compareResults.Fisheries.SelFcns(models[compare],fleets="TCF",singlePlot=TRUE,
                                            dodge=0,years=c(1990,2013:2024))
rCompTCMs::compareResults.Fisheries.SelFcns(models[compare],fleets="GF All",singlePlot=TRUE,
                                            dodge=0,years=c(1985,1995,2015))

rCompTCMs::compareFits.MaturityOgiveData(models[compare],
                                         dodge=0,
                                         types="fits",
                                         plot1stObs=FALSE);
#
# rCompTCMs::compareFits.EffectiveNs(models[compare],
#                                    fleet.type="survey",
#                                    facet_grid="x+m+s~.",
#                                    dodge=0)
#
# rCompTCMs::compareFits.EffectiveNs(models[compare],
#                                    fleet.type="fishery",
#                                    category="retained",
#                                    facet_grid="x+m+s~.",
#                                    dodge=0)
#
# rCompTCMs::compareFits.EffectiveNs(models[compare],
#                                    fleet.type="fishery",
#                                    category="total",
#                                    facet_grid="x+m+s~.",
#                                    dodge=0)
#
if (FALSE){
  #--plot fits to abundance, biomass data
  #----surveys
  rCompTCMs::compareFits.BiomassData(models[compare],fleets="NMFS M",
                                     fleet.type="survey",catch.type="index",
                                     plot1stObs=FALSE,numRecent=30)
  rCompTCMs::compareFits.BiomassData(models[compare],fleets="NMFS F",
                                     fleet.type="survey",catch.type="index",
                                     plot1stObs=FALSE,numRecent=30)
  rCompTCMs::compareFits.AbundanceData(models[compare],fleet.type="survey",catch.type="index",
                                     plot1stObs=FALSE,numRecent=30)
  #----fisheries
  rCompTCMs::compareFits.BiomassData(models[compare],fleet.type="fishery",catch.type="retained",
                                     fishery.pdfType = "lognormal",
                                     plot1stObs=FALSE,numRecent=20)
  rCompTCMs::compareFits.BiomassData(models[compare],fleet.type="fishery",catch.type="total",
                                     fishery.pdfType = "lognormal",
                                     plot1stObs=FALSE,numRecent=20)
  rCompTCMs::compareFits.AbundanceData(models[compare],fleet.type="fishery",catch.type="retained",
                                     fishery.pdfType = "lognormal",
                                     plot1stObs=FALSE,numRecent=20)
  rCompTCMs::compareFits.AbundanceData(models[compare],fleet.type="fishery",catch.type="total",
                                     fishery.pdfType = "lognormal",
                                     plot1stObs=FALSE,numRecent=20)
  
  #--plot fits to size composition data
  #----surveys
  rCompTCMs::compareFits.MeanSizeComps(models[compare],fleets="NMFS M",
                                       fleet.type="survey",catch.type="index",
                                       plot1stObs=TRUE);
  rCompTCMs::compareFits.MeanSizeComps(models[compare],fleets="NMFS F",
                                       fleet.type="survey",catch.type="index",
                                       plot1stObs=TRUE);
  compareFits.SizeComps(models[compare],fleets="NMFS M",
                                       fleet.type="survey",catch.type="index",
                                       plot1stObs=FALSE);
  compareFits.SizeComps(models[compare],fleets="NMFS F",
                                       fleet.type="survey",catch.type="index",
                                       plot1stObs=FALSE);
  #----fisheries
  rCompTCMs::compareFits.MeanSizeComps(models[compare],
                                       fleet.type="fishery",catch.type="total",
                                       plot1stObs=TRUE);
  compareFits.SizeComps(models[compare],fleets="TCF",
                                       fleet.type="fishery",catch.type="retained",
                                       plot1stObs=TRUE);
  compareFits.SizeComps(models[compare],fleets="TCF",
                                       fleet.type="fishery",catch.type="total",
                                       plot1stObs=TRUE);
  compareFits.SizeComps(models[compare],fleets="SCF",
                                       fleet.type="fishery",catch.type="total",
                                       plot1stObs=TRUE);
  compareFits.SizeComps(models[compare],fleets="GF All",
                                       fleet.type="fishery",catch.type="total",
                                       plot1stObs=TRUE);
  
  #--plot fits to other data
  rCompTCMs::compareFits.MaturityOgiveData(models[compare],
                                           dodge=0,
                                           types="fits",
                                           plot1stObs=TRUE);
  rCompTCMs::compareFits.GrowthData(models[compare],
                                           dodge=0,
                                           plot1stObs=TRUE);
  
  #--plot population results
  library(ggplot2);
  dfrNMs<-rTCSAM02::getMDFR.Pop.NaturalMortality(models[compare],type="M_cxm");
  dfrNMsW.Imm<-reshape2::dcast(dfrNMs[dfrNMs$m=="IMMATURE",],y+case~x,fun.aggregate=wtsUtilities::Sum,value.var="val");
  dfrNMsW.Mat<-reshape2::dcast(dfrNMs[dfrNMs$m=="MATURE",],  y+case~x,fun.aggregate=wtsUtilities::Sum,value.var="val");
  rCompTCMs::compareResults.Pop.NaturalMortality.BarPlot(models[compare]);
  rCompTCMs::compareResults.Pop.NaturalMortality(models[compare],dodge=0)
  rCompTCMs::compareResults.Pop.MeanGrowth(models[compare],dodge=0)
  rCompTCMs::compareResults.Pop.GrowthMatrices.LinePlots(models[compare],dodge=0)
  rCompTCMs::compareResults.Pop.PrM2M(models[compare],dodge=0)
  rCompTCMs::compareResults.Pop.Recruitment(models[compare],dodge=0,numRecent=30)
  rCompTCMs::compareResults.Pop.MatureBiomass(models[compare],dodge=0,numRecent=30)
  rCompTCMs::compareResults.Pop.Abundance(models[compare],cast="x+m",dodge=0,years="all",facet_grid=x+m~.);
  rCompTCMs::compareResults.Pop.CohortProgression(models[compare],years=1:10,cast="x+m+z",
                                                  facet_grid=x+m~.,scales="free_y");
  
  #--plot surveys-related quantities
  rCompTCMs::compareResults.Surveys.Catchability(models[compare],dodge=0)
  rCompTCMs::compareResults.Surveys.SelFcns(models[compare],singlePlot=TRUE,
                                            dodge=0,years=c(1981,2018));
  rCompTCMs::compareResults.Surveys.CaptureProbs(models[compare],singlePlot=TRUE,
                                                 dodge=0,years=c(1981,2018));
  
  #--plot fisheries-related quantities
  rCompTCMs::compareResults.Fisheries.Catchability(models[compare],dodge=0)
  rCompTCMs::compareResults.Fisheries.RetFcns(models[compare],fleets="TCF",singlePlot=TRUE,
                                              dodge=0,years=c(1980,1991,2005,2007,2009,2014,2017))
  rCompTCMs::compareResults.Fisheries.SelFcns(models[compare],fleets="TCF",singlePlot=TRUE,
                                              dodge=0,years=c(1980,1991:2018))
  rCompTCMs::compareResults.Fisheries.SelFcns(models[compare],fleets="SCF",singlePlot=TRUE,
                                              dodge=0,years=c(1990,2000,2015))
  rCompTCMs::compareResults.Fisheries.SelFcns(models[compare],fleets="RKF",singlePlot=TRUE,
                                              dodge=0,years=c(1990,2000,2015))
  rCompTCMs::compareResults.Fisheries.SelFcns(models[compare],fleets="GF All",singlePlot=TRUE,
                                              dodge=0,years=c(1980,1990,2000))
}
#
#create pdf output
rCompTCMs::modelComparisons.ModelFits.OtherData(models[compare],
                                               plot1stObs=TRUE,
                                               output_format="pdf_document",
                                               clean=TRUE);
rCompTCMs::modelComparisons.ModelFits.ACD.Fisheries(models[compare],
                                         plot1stObs=FALSE,
                                         output_format="pdf_document",
                                         clean=TRUE);
# rCompTCMs::modelComparisons.ModelFits.ACD.Surveys(models[compare],
#                                          plot1stObs=TRUE,
#                                          output_format="pdf_document",
#                                          clean=TRUE);
# rCompTCMs::modelComparisons.ModelFits.ZCsByYear(models[compare],
#                                          plot1stObs=TRUE,
#                                          output_format="pdf_document",
#                                          clean=TRUE);
# rCompTCMs::modelComparisons.ModelFits.ZCs(models[compare],
#                                           type="Surveys",
#                                           plot1stObs=TRUE,
#                                           output_format="pdf_document",
#                                           clean=TRUE);
# rCompTCMs::modelComparisons.ModelFits.ZCs(models[compare],
#                                           type="Fisheries",
#                                           plot1stObs=TRUE,
#                                           output_format="pdf_document",
#                                           clean=TRUE);
# rCompTCMs::modelComparisons.PopProcesses(models[compare],
#                                          output_format="pdf_document",
#                                          clean=TRUE);
rCompTCMs::modelComparisons.PopQuantities(models[compare],
                                         output_format="pdf_document",
                                         clean=TRUE);
# ##--NMFS EBS
# rCompTCMs::modelComparisons.Characteristics.Surveys(models[compare],
#                                                     surveys=c("NMFS M","NMFS F"),
#                                                     selYears=c(1975, 1988),
#                                                     avlYears=NULL,
#                                                     capYears=c(1975,1988),
#                                                     output_format="pdf_document",
#                                                     clean=TRUE);
# ##--SBS
# rCompTCMs::modelComparisons.Characteristics.Surveys(models[compare],
#                                                     surveys=c("SBS BSFRF males",
#                                                               "SBS BSFRF females",
#                                                               "SBS NMFS males",
#                                                               "SBS NMFS females"),
#                                                     selYears=c(2013),
#                                                     avlYears=c(2013,2014,2015,2016,2017),
#                                                     capYears=c(2013,2014,2015,2016,2017),
#                                                     output_format="pdf_document",
#                                                     clean=TRUE);
# ##--TCF (need to rename output pdf)
# rCompTCMs::modelComparisons.Characteristics.Fisheries(models[compare],
#                                                        fisheries="TCF",
#                                                        selyears=c(1990,1991:2018),
#                                                        retyears=c(1990,1995,2018),
#                                                        output_format="pdf_document",
#                                                        clean=TRUE);
# ##--SCF
# rCompTCMs::modelComparisons.Characteristics.Fisheries(models[compare],
#                                                        fisheries="SCF",
#                                                        selyears=c(1995,2000,2015),
#                                                        retyears=NULL,
#                                                        output_format="pdf_document",
#                                                        clean=TRUE);
# ##--RKF
# rCompTCMs::modelComparisons.Characteristics.Fisheries(models[compare],
#                                                        fisheries="RKF",
#                                                        selyears=c(1990,2000,2015),
#                                                        retyears=NULL,
#                                                        output_format="pdf_document",
#                                                        clean=TRUE);
# ##--GTF
# rCompTCMs::modelComparisons.Characteristics.Fisheries(models[compare],
#                                                        fisheries="GTF",
#                                                        selyears=c(1985,1995,2015),
#                                                        retyears=NULL,
#                                                        output_format="pdf_document",
#                                                        clean=TRUE);
# rCompTCMs::modelComparisons.ParameterTables(models[compare],
#                                            output_format="pdf_document");

# rCompTCMs::modelComparisons(models[compare],
#                             plot1stObs=FALSE,
#                             output_format="pdf_document");
library(ggplot2);
dfrNMs<-rTCSAM02::getMDFR.Pop.NaturalMortality(models[compare],type="M_cxm");
dfrNMsW.Imm<-reshape2::dcast(dfrNMs[dfrNMs$m=="IMMATURE",],y+case~x,fun.aggregate=wtsUtilities::Sum,value.var="val");
dfrNMsW.Mat<-reshape2::dcast(dfrNMs[dfrNMs$m=="MATURE",],  y+case~x,fun.aggregate=wtsUtilities::Sum,value.var="val");
mxNM<-max(dfrNMs$val);
uMs<-unique(tolower(dfrNMs$m));
plots<-list();
for (uM in uMs){
  tmp<-dfrNMs[tolower(dfrNMs$m)==uM,];
  p <- ggplot(tmp,aes_string(x="x",y="val",fill="case"))+geom_bar(stat="identity",position=position_dodge());
  p <- p + facet_grid(y~m);
  p <- p + ylim(0,mxNM);
  p <- p + labs(x="",y="natural mortality",fill="scenario");
  pl<- p; #--save plot w/ legend
  p <- p + theme(legend.position="none",);
  plots[[uM]]<-p;
}
prows<-cowplot::plot_grid(plotlist=plots,nrow=2,rel_heights=c(1,2));
legend<-cowplot::get_legend(pl);
p<-cowplot::plot_grid(prows,legend,ncol=2,rel_widths=c(4,2));
pdf("NMs.pdf",width=8,height=6);
print(p);
dev.off();

