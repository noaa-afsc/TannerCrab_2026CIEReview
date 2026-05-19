#convert TCSAM02 model output to RData files
require(rTCSAM02)

dirThis = dirname((rstudioapi::getSourceEditorContext())$path);
setwd(dirThis);

case<-c("22_03d5d");
runDir = "run";
resDir = "./run_results";
if (!dir.exists(resDir)) dir.create(resDir);

#--remove old results files
fns1<-list.files(path=resDir,pattern=glob2rx("OFCs*.csv"));
fns2<-list.files(path=resDir,pattern=glob2rx("Results*.RData"));
if (length(fns1)>0) file.remove(file.path(resDir,fns1))
if (length(fns2)>0) file.remove(file.path(resDir,fns2));

#--create new results files
best<-rTCSAM02::getResLst(file.path(runDir),verbose=TRUE)
if (!is.null(best)) {
  wtsUtilities::saveObj(best,fn=file.path(resDir,"Results.RData"));
  mdfr<-rTCSAM02::getMDFR.OFCs.DataComponents(best,verbose=FALSE);
  readr::write_csv(mdfr,file=file.path(resDir,"OFCs.DataComponents.csv"));
  mdfrp<-rTCSAM02::getMDFR.OFCs.NonDataComponents(best,verbose=FALSE)
  readr::write_csv(mdfrp,file=file.path(resDir,"OFCs.NonDataComponents.csv"));
  rm(best,mdfr,mdfrp)
}

models<-list();
models[[case]] = wtsUtilities::getObj(fn=file.path(resDir,"Results.RData"));
compare = names(models);

#--extract tables and create plots
dfr.PVs   <-rCompTCMs::extractMDFR.Results.ParameterValues(models[compare]);
readr::write_csv(dfr.PVs,file.path(resDir,"ParamValues.csv"))
dfr.PsAtBs<-rCompTCMs::extractMDFR.Results.ParametersAtBounds(models[compare],delta=0.001);
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

if (file.exists(file.path(runDir,"admodel.hes"))){
  adnuts::check_identifiable("tcsam02",path=runDir)
  lst = adnuts:::.read_mle_fit("tcsam02",path=runDir)
  hes = adnuts:::.getADMBHessian(path=runDir)
  evs_hes = eigen(hes)$values
  min_hes = head(sort(evs_hes)); print(min_hes);
}


