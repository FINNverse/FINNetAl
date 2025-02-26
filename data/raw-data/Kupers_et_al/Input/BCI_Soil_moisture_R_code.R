###.......................................................................###
### This code accompanies the Data Descriptor by Kupers et al. titled     ###
### 'Dry season soil water potential maps of a 50 hectare tropical forest ###
### plot on Barro Colorado Island, Panama' (2019) Scientific Data.        ###
###   N.B. Input data required for running the code is available through  ###
###   the links provided in the section 'Load packages and data' below.   ###
###.......................................................................###

#.......................................................#
#### Set working directories and map settings        ####
#.......................................................#

###Set folder containing R code (current folder) as input
wd.input = dirname(rstudioapi::getActiveDocumentContext()$path)
###Set output folder
wd.output = paste0(unlist(strsplit(wd.input, split='Input', fixed=TRUE)),"Output")

###Run the Random Forest model?
#If FALSE, the model provided in the .RData object will be loaded.
run.model = F

###Run custom map and/or the standard maps?
custom.map = T
standard.maps = F
###Map settings for custom map
#Set date or monitored soil water content (SWC)
use.date = T
if(use.date == T){set.date = "2016-03-01"}  #Set any date starting from the year 1975. Use the format "YYYY-MM-DD".
if(use.date == F){set.swc = 45}             #Monitored soil water content, i.e. dry season stage. Median of measurements: XX %.

#Other settings for modeling and mapping
set.time = 0.5                              #Set time for custom maps. Fraction of 24 hours, e.g. 0.375 = 9 AM, 0.5 = 12 PM. Median of measurements:: 0.47.
set.depth = 15                              #Set depth for custom maps. Median of measurements: 16 cm.
set.ba = 0.03                               #Set basal area for custom maps. Median of the 2015 census: 0.03 m 25m^-2
tuned = F                                   #When running the Random Forest model, tune the nodesize and mtry parameters or not?
if(tuned==F){nodesize = 5; mtry = 3}        #defaults are nodesize = 5; mtry = 3
ntree = 1000                                #default is 1000

#Save plots and tables as presented in the Data Descriptor?
save.figures = F


#.......................................................#
#### Load packages and data                          ####
#.......................................................#

###Packages
library(lubridate)
library(randomForestSRC)
library(raster)
library(remotes) #or: library(devtools)
library(bciex)
library(RSAGA)
library(caret)
library(zoo)
library(vegan)
library(ModelMetrics)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(ggsci)
library(plotrix)
library(xlsx)

###Data and additional functions
setwd(wd.input)
load("CTFSRPackage.RData")                                                #Available via: http://ctfs.si.edu/Public/CTFSRPackage/
moisture = read.table("BCI_Soil_moisture_mapping.txt", header = T)        #Available via: Data citation 1: Kupers, S.J. et al. Dry season soil water potential maps of a 50 hectare tropical forest plot on Barro Colorado Island, Panama. Scientific Data.
small.scale = read.table("BCI_Soil_moisture_small_scale.txt",header=T)    #Available via: Data citation 1: Kupers, S.J. et al. Dry season soil water potential maps of a 50 hectare tropical forest plot on Barro Colorado Island, Panama. Scientific Data.
monitored.swc.stri = read.csv("bci_lutz_shallow_gsm_man.csv")             #Available via: https://biogeodb.stri.si.edu/physical_monitoring/research/barrocolorado
monitored.rain.stri = read.csv("bci_cl_ra_man.csv")                       #Available via: https://biogeodb.stri.si.edu/physical_monitoring/research/barrocolorado
topography = bci_elevation                                                #Available via:  bciex package (https://rdrr.io/github/forestgeo/bciex/)
habitats = bci_habitat                                                    #Available via:  bciex package (https://rdrr.io/github/forestgeo/bciex/). Original citation: Harms et al. (2001). Habitat associations of trees and shrubs in a 50‐ha Neotropical forest plot. Figure 1. Journal of Ecology 89, 947-959, doi:https://doi.org/10.1111/j.1365-2745.2001.00615.x
soiltype = read.table("BCI_soil_type.txt", sep = '\t',header=T)           #Available via: Data citation 1: Kupers, S.J. et al. Dry season soil water potential maps of a 50 hectare tropical forest plot on Barro Colorado Island, Panama. Scientific Data. Original citation: Baillie et al. (2007). Semi-detailed soil survey of Barro Colorado Island, Panama. Figure 7.1. https://biogeodb.stri.si.edu/bioinformatics/bci_soil_map/
if(run.model==T){
  census2015 = read.table("PlotDataReport03-28-2018_890930496.txt",       #Available via: http://ctfs.si.edu/webatlas/datasets/bci/
                          sep = '\t',header=T)        
}
#Load the Random Forest model object when only mapping 
if(run.model == F){load("BCI_SWP_RF_model.RData")}                        #Available via: Data citation 1: Kupers, S.J. Dry season soil water potential maps of the 50-ha Forest Dynamics Plot on BCI, Panama. iDiv Data Repository.


#.......................................................#
#### Prepare input data                              ####
#.......................................................#


###Elevation
elevation = topography[,c("x","y","elev")]

###Slope
elev.r = rasterFromXYZ(elevation)
projection(elev.r) ##The raster needs a projection
crs(elev.r)  <- "+proj=lcc +lat_1=9.15+lon_0=-79.85 +ellps=WGS84"
slope.r = terrain(elev.r, opt='slope', unit='degrees', neighbors=8)
slope = as.data.frame(slope.r,xy=T)
slope$slope = slope.r@data@values
slope = slope[complete.cases(slope),]

###Soil type
soiltype$soiltype = as.numeric(soiltype$soiltype)  #1 = AVA, 2 = Fairchild, 3 = Marron, 4 = Swamp
##Adjust the swamp location with the habitat maps from Harms et al. (2001)
colnames(habitats)[3] = "habitat.name"
habitats.tmp = data.frame(habitat.name = unique(habitats$habitat.name),habitat = 1:length(unique(habitats$habitat.name)))
habitats = merge(habitats,habitats.tmp)
habitats = habitats[,-which(colnames(habitats)=="habitat.name")]
habitats.r = rasterFromXYZ(habitats) #swamp = habitat 5
soiltype.tmp = pick.from.points(data=soiltype, src = habitats, method = c("nearest.neighbour"), X.name = "x", Y.name = "y", cbind=TRUE)
soiltype.tmp$soiltype = ifelse(soiltype.tmp$soiltype==4&soiltype.tmp$habitat!=5,3,soiltype.tmp$soiltype)
soiltype.tmp$soiltype = ifelse(soiltype.tmp$habitat==5,4,soiltype.tmp$soiltype)
soiltype.tmp = soiltype.tmp[,c("x","y","soiltype")]
soiltype.r = rasterFromXYZ(soiltype.tmp)
soiltype = soiltype.tmp

###Monitored soil water content
monitored.swc = monitored.swc.stri[which(monitored.swc.stri$chk_note=="good"),]
monitored.swc = aggregate(monitored.swc$h2o_by_dry,list(monitored.swc$date),mean,na.rm=T)
colnames(monitored.swc) = c("date","swc")
monitored.swc$date = format(as.POSIXct(strptime(monitored.swc$date,"%d/%m/%Y",tz="")) ,format = "%Y-%m-%d")
monitored.swc = monitored.swc[order(monitored.swc$date),]
##Add all dates
all.dates = data.frame(date = seq(as.Date(monitored.swc$date[1]), as.Date(monitored.swc$date[length(monitored.swc$date)]), by="days"))
all.dates$date = as.character(all.dates$date)
monitored.swc$date = as.character(monitored.swc$date)
monitored.swc = merge(all.dates,monitored.swc,all.x=T)
##Linearly interpolate values for missing dates
swc.filled = zoo(monitored.swc$swc)
swc.filled = na.approx(swc.filled)
swc.filled = as.numeric(swc.filled)
monitored.swc$swc.filled = swc.filled
monitored.swc.days = unique(moisture$date)
monitored.swc.days = format(as.POSIXct(strptime(monitored.swc.days,"%m/%d/%Y",tz="")) ,format = "%Y-%m-%d")
monitored.swc.days = as.character(monitored.swc.days)
monitored.swc.days = data.frame(date = monitored.swc.days)
monitored.swc.days = merge(monitored.swc.days,monitored.swc[,c("date","swc.filled")],all.x=T)
#Calculate monitored SWC in an early, mid and late dry season (25th, 50th and 75th percentile SWC in February until April across years.
swc.regular = monitored.swc[,c("date","swc.filled")]
swc.regular$month = format(as.POSIXct(strptime(swc.regular$date,"%Y-%m-%d",tz="")) ,format = "%m")
swc.regular$month = as.numeric(as.character(swc.regular$month))
swc.regular$day = format(as.POSIXct(strptime(swc.regular$date,"%Y-%m-%d",tz="")) ,format = "%d")
swc.regular$day = as.numeric(as.character(swc.regular$day))
swc.regular = swc.regular[which(swc.regular$month%in%2:4),]
swc.regular = quantile(swc.regular$swc.filled, probs = seq(0.25,0.75,length.out=3),na.rm=T)
#Calculate monitored SWC during a drought (the median of SWC in the March 2016 measurement period)
swc.elnino = monitored.swc.days[grepl("2016",monitored.swc.days$date),]
swc.elnino = median(swc.elnino$swc.filled)
#Merge for later map making
swc.periods = data.frame(period = c("Early (regular year)","Progressing (regular year)","Peak (regular year)","Peak (drought year)"),
                         monitored.swc = c(rev(swc.regular),swc.elnino))

###Basal area
if(run.model==T){
  census2015.ba = census2015[which(census2015$Status=="alive"),]
  census2015.ba$ba = ba(census2015.ba$DBH,dbhunit = "mm")
  census2015.ba = census2015.ba[,c("PX","PY","ba")]
  census2015.ba = census2015.ba[complete.cases(census2015.ba),]
  census2015.ba.r = rasterize(x = census2015.ba$PX, y = census2015.ba$PY, z = census2015.ba$ba,
                              gridsize=5, plotdim=c(1000,500),graph=F)
  coordinates = expand.grid(seq(2.5,497.5,5),seq(2.5,997.5,5))
  colnames(coordinates) = c("y","x")
  census2015.ba = c(census2015.ba.r)
  census2015.ba = cbind(coordinates,census2015.ba)
  census2015.ba$census2015.ba = log(census2015.ba$census2015.ba)
  
}

###Merge all predictors with moisture data and select input variables
if(run.model==T){
  input = moisture
  input$depth = log(input$depth)
  input$date = format(as.POSIXct(strptime(input$date,"%m/%d/%Y",tz="")) ,format = "%Y-%m-%d")
  input = merge(input,monitored.swc.days)
  input = pick.from.points(data=input, src = elevation, method = c("nearest.neighbour"),X.name = "x", Y.name = "y", cbind=TRUE)
  input = pick.from.points(data=input, src = slope, method = c("nearest.neighbour"), X.name = "x", Y.name = "y", cbind=TRUE)
  input = pick.from.points(data=input, src = soiltype, method = c("nearest.neighbour"), X.name = "x", Y.name = "y", cbind=TRUE)
  input$soiltype = as.factor(input$soiltype)
  input = pick.from.points(data=input, src = census2015.ba, method = c("nearest.neighbour"),X.name = "x", Y.name = "y", cbind=TRUE)

  input.period = input[,c("period","swp","col.time","swc.filled","elev","slope","soiltype","census2015.ba","depth")] # to check performance per period.
  input = input[,c("swp","col.time","swc.filled","elev","slope","soiltype","census2015.ba","depth")] # for running the RF model.
  pred.xlabels = c("Time (hh:mm)","Monitored SWC (%)","Elevation (m)","Slope (°)","Soil type",as.expression(bquote("Basal area (m"^2*" 25 m"^-2*")")),"Depth (cm)")
}


#.......................................................#
#### Run and validate Random Forest model            ####
#.......................................................#


#Run model and plot performance and predictors if wanted (set in line 19)
setwd(wd.output)
if(run.model==T){
  ###Run model
  #Set in line 35 if tuning of the model settings should be done.
  #We use default settings below
  if(tuned==1){
    rf.tune = tune(formula = swp ~ ., data = input)
    nodesize = rf.tune$optimal[1]; mtry = rf.tune$optimal[2]
  }else{
    nodesize = nodesize; mtry = mtry  #Set in line 34
  }
  ##Best settings of tuned model used
  rf.model = rfsrc(formula = swp ~ ., data = input, nodesize = nodesize, mtry = mtry, 
                   importance=T,na.action = "na.impute",seed=-1,ntree = ntree)
  rf.model
  save(rf.model, file = paste("BCI_SWP_RF_model_",Sys.Date(),".RData",sep="")) #Save with date, as to not overwrite other models. If needed replace Sys.date() with Sys.time().

  
  ###Evaluate model performance
  if(save.figures==T){
    pdf("Figure 5 - Model performance.pdf",width=12,height=10)
  }
  par(mfrow = c(1,2))
  layout(rbind(c(rep(1,2),rep(6,2)),c(rep(1,2),rep(6,2)), c(2,3,7,8),c(4,5,9,10)))
  periods = unique(input.period$period)
  periods.names.left = c("c) Feb 2015","d) Mar 2015","g) Apr 2015","h) Mar 2016")
  periods.names.right = c("e) Feb 2015","f) Mar 2015","i) Apr 2015","j) Mar 2016")
  colors = pal_npg(palette = c("nrc"), alpha = 0.4)(4)
  #Plot
  for(i in 1:2){
    for(j in 1:5){
      if(i==1){
        if(j==1){
          obs = input$swp
          pred = rf.model$predicted.oob
          mtext = "a) Out-of-bag performance"
          xlab = "Predicted SWP (MPa)"
          ylab = "Observed SWP (MPa)"
          color = rgb(.2,.2,.2,alpha=.4)
          par(mar=c(5,5,4,1))
          cex.factor = 1
          leg.dist = c(0.05,0.11,0.17)
          plot.xlab = T; plot.ylab = T
        }else{
          period = periods[j-1]
          rows = which(input.period$period==period)
          obs = input.period[rows,"swp"]
          pred = rf.model$predicted.oob[rows]
          mtext = periods.names.left[j-1]
          cex.factor = 0.8
          color = colors[j-1]
          leg.dist = c(0.05,0.15,0.23)
          if(j==2){par(mar=c(1,5,6,0)); plot.xlab = F; plot.ylab = T; xlab = NA; ylab = "Observed SWP (MPa)"}
          if(j==3){par(mar=c(1,4,6,1)); plot.xlab = F; plot.ylab = F; xlab = NA; ylab = NA}
          if(j==4){par(mar=c(5,5,2,0)); plot.xlab = T; plot.ylab = T; xlab = "Predicted SWP (MPa)"; ylab = "Observed SWP (MPa)"}
          if(j==5){par(mar=c(5,4,2,1)); plot.xlab = T; plot.ylab = F; xlab = "Predicted SWP (MPa)"; ylab = NA}
        }
      }else{
        if(j==1){
          obs = input$swp
          pred = rf.model$predicted
          mtext = "b) Full model predictions"
          xlab = "Predicted SWP (MPa)"
          ylab = NA
          color = rgb(.2,.2,.2,alpha=.4)
          par(mar=c(5,5,4,1))
          cex.factor = 1
          leg.dist = c(0.05,0.11,0.17)
          plot.xlab = T; plot.ylab = T
        }else{
          period = periods[j-1]
          rows = which(input.period$period==period)
          obs = input.period[rows,"swp"]
          pred = rf.model$predicted[rows]
          mtext = periods.names.right[j-1]
          cex.factor = 0.8
          color = colors[j-1]
          leg.dist = c(0.05,0.15,0.23)
          if(j==2){par(mar=c(1,5,6,0)); plot.xlab = F; plot.ylab = T; xlab = NA; ylab = NA}
          if(j==3){par(mar=c(1,4,6,1)); plot.xlab = F; plot.ylab = F; xlab = NA; ylab = NA}
          if(j==4){par(mar=c(5,5,2,0)); plot.xlab = T; plot.ylab = T; xlab = "Predicted SWP (MPa)"; ylab = NA}
          if(j==5){par(mar=c(5,4,2,1)); plot.xlab = T; plot.ylab = F; xlab = "Predicted SWP (MPa)"; ylab = NA}
        }
      }
      ##Performance
      r2 = 1 - sum((pred - obs) ^ 2)/sum((obs - mean(obs)) ^ 2)
      rmse = rmse(obs,pred)
      mae = mae(obs,pred)
      #rmse.mae = rmse/mae
      r2 = sprintf("%.2f", r2); r2
      rmse = sprintf("%.2f", rmse); rmse
      mae = sprintf("%.2f", mae); mae
      
      #Plot points
      xlim = c(-1.75,0)
      ylim = c(-2.5,0)
      plot(obs~pred,xlim = xlim, ylim = ylim, xaxt="n", yaxt="n",xlab = xlab, ylab = ylab,
           pch=19, col = "white", cex.lab=1.5*cex.factor)  #,cex.axis=1.5*cex.factor
      if(plot.xlab == T){axis(side=1,cex.axis=1.5*cex.factor)}else{axis(side=1,labels=NA)}
      if(plot.ylab == T){axis(side=2,cex.axis=1.5*cex.factor)}else{axis(side=2,labels=NA)}
      mtext(mtext, side=3, line=0.5, at=-1.75,cex=2*cex.factor,adj=0)
      abline(0,0,col="lightgrey");abline(v=0,col="lightgrey")
      abline(0,1,lty=2)
      points(obs~pred,pch=19,col=color)
      text(-1.75,-1.85, labels = "1:1")
      #Plot values
      label = c(as.expression(bquote("R"^2 == .(r2))), bquote("RMSE" == .(rmse)), bquote("MAE" == .(mae)))
      text(rep(xlim[1]+(0.00*diff(xlim)),length(label)),
           c(ylim[2]-(leg.dist[1]*diff(ylim)),ylim[2]-(leg.dist[2]*diff(ylim)),ylim[2]-(leg.dist[3]*diff(ylim))),
           labels=label,adj=0,cex=1.4*cex.factor)
      #Plot regression line
      clip(x1=range(pred)[1],x2=range(pred)[2],y1=range(obs)[1], y2=range(obs)[2])
      abline(lm(obs~pred),lwd=2)
    }
  }
  if(save.figures==T){dev.off()}
  
  ###Plot fitted SWP versus predictors of the Random Forest model
  if(save.figures==T){
    pdf("Figure 6 - Fitted SWP vs predictors.pdf",width = 10,height = 5)
  }
  par(mfrow = c(2,4))
  fitted.input = input[,-1]
  ##First set all predictors to their mean
  #Determine the most common soil type (as mean value when plotting other predictors)
  common.soil = table(fitted.input$soiltype)
  common.soil = names(common.soil[which.max(common.soil)])
  fitted.input$soiltype = as.numeric(as.character(fitted.input$soiltype))
  #Get means of other values
  fitted.input = colMeans(fitted.input,na.rm=T)
  fitted.input[which(names(fitted.input)=="soiltype")] <- as.numeric(common.soil)
  fitted.input = data.frame(as.list(fitted.input))
  fitted.input = fitted.input[rep(seq_len(nrow(fitted.input)), each=100),]
  fitted.input = as.data.frame(lapply(fitted.input, as.numeric))
  ##Plot
  for(i in 1:length(fitted.input)){
    fitted.input.tmp = fitted.input
    #Determine importance and scale to 100%
    imp = rf.model$importance
    imp = imp / sum(imp) * 100
    imp = imp[order(imp,decreasing=T)]
    ##Set focal variable to its range
    xvar.col = which(colnames(fitted.input.tmp)==names(imp)[i])
    xvar.tmp = colnames(fitted.input.tmp)[xvar.col]
    xlab.tmp = pred.xlabels[xvar.col]
    if(xvar.tmp!="soiltype"){
      fitted.input.tmp[,xvar.tmp] = seq(min(input[,xvar.tmp],na.rm=T),max(input[,xvar.tmp],na.rm=T),length=100)
      fitted.input.tmp$soiltype = as.factor(fitted.input.tmp$soiltype)
    }else{
      fitted.input.tmp = fitted.input.tmp[!duplicated(fitted.input.tmp),] 
      fitted.input.tmp = fitted.input.tmp[rep(seq_len(nrow(fitted.input.tmp)), each=length(unique(input$soiltype))),]
      fitted.input.tmp$soiltype = as.factor(seq(1,length(unique(input$soiltype)),1))
    }
    #Predict SWP and plot
    pred = predict(rf.model, fitted.input.tmp)$predicted
    plot(pred~fitted.input.tmp[,xvar.tmp],type="l",ylab="Predicted SWP (MPa)",xaxt="n",xlab=xlab.tmp,cex.lab=1.2)
    #Plot labels
    if(!xvar.tmp%in%c("soiltype","col.time","depth","census2015.ba")){axis(side = 1,at=NULL)}
    if(xvar.tmp=="soiltype"){axis(side = 1,at=c(1,2,3,4),labels = c("A","F","M","S"))}
    #if(xvar.tmp=="soiltype"){axis(side = 1,at=c(1,2,3,4,5),labels = c("Av","Fa","Ma","Sw","Mx"))}
    if(xvar.tmp=="col.time"){axis(side = 1,at=c(0.375,0.500,0.625),labels = c("09:00","12:00","15:00"))}
    if(xvar.tmp=="depth"){axis(side = 1,at=log(c(1,10,25,50,100)),labels=c(1,10,25,50,100))}
    if(xvar.tmp=="census2015.ba"){axis(side = 1,at=log(c(0.05,0.2,0.5,1,2,4)),labels=c(0.05,0.2,0.5,1,2,4))}
    ##Add quantiles
    if(xvar.tmp!="soiltype"){
      quantiles = quantile(input[,xvar.tmp], probs = seq(0,1,0.1),na.rm=T)
      segments(x0 = quantiles, x1 = quantiles,
               y0 = max(pred) + 0.05*diff(range(pred)), y1 = max(pred) + 0.1*diff(range(pred)),xpd = TRUE)
    }
    ##Add importance value
    imp.round = sprintf("%.2f", imp[i])
    mtext(paste0(letters[i],") Importance = ",imp.round,"%"), side=3, line=1.0,
          at=par("usr")[1]+0.05*diff(par("usr")[1:2]),cex=0.80,adj=0)
  }
  if(save.figures==T){dev.off()}
}



#.......................................................#
#### Construct soil moisture maps                    ####
#.......................................................#

###Set SWP scale and colors
zlim = c(-1.30,0)
colors = colorRampPalette(c('darkred','red','orange','yellow','lightblue','blue'))(50)

###Prepare input data on grid
x.pred = seq(2.5,997.5,5)
y.pred = seq(2.5,497.5,5)
pred.map = expand.grid(x.pred,y.pred)
colnames(pred.map) = c("x","y")
#Fixed values
pred.map$col.time = set.time; set.time
pred.map$depth = set.depth; set.depth
pred.map$census2015.ba = set.ba; set.ba
#Set monitored soil water content, optionally for a specific date (set in line 26-27)
if(use.date==T){set.swc = monitored.swc[which(monitored.swc$date==set.date),"swc.filled"]}
pred.map$swc.filled = set.swc;set.swc
#Spatial variables
pred.map = pick.from.points(data=pred.map, src = elevation, method = c("nearest.neighbour"), X.name = "x", Y.name = "y", cbind=TRUE)
pred.map = pick.from.points(data=pred.map, src = slope, method = c("nearest.neighbour"), X.name = "x", Y.name = "y", cbind=TRUE)
pred.map = pick.from.points(data=pred.map, src = soiltype, method = c("nearest.neighbour"),X.name = "x", Y.name = "y", cbind=TRUE)
pred.map$soiltype = as.factor(pred.map$soiltype)

###Plot and save custom map, as set in line 22
if(custom.map==T){
  if(use.date==T){
    filename = paste0("BCI_SWP_map_date_=_",set.date)
  }else{
    filename = paste0("BCI_SWP_map_monitored_SWC_=_",round(set.swc,2),"_percent")
  }
  swp = predict(rf.model, pred.map)$predicted
  pred.swp = cbind(pred.map[,c("x","y")],swp)
  pred.swp$swp = ifelse(pred.swp$swp>0,0,pred.swp$swp)
  #.txt
  write.table(pred.swp,paste(filename,"txt",sep="."), row.names=F,sep="\t")
  #.pdf
  swp.r = rasterFromXYZ(pred.swp)
  pdf(paste(filename,".pdf",sep=""),width = 10, height = 5.69)
  plot(swp.r, col = colors,zlim = zlim,
       smallplot=c(0.88,0.9, 0.3,0.75),
       legend.args=list(text='MPa', side=1, font=1, line=-15.0, cex=1,adj=0))
  elevation.m = matrix(pred.map$elev,ncol = length(x.pred), nrow = length(y.pred),byrow=T)
  maptopo(elevation.m,plotdim=c(1000,500),add=T, interval=2)
  dev.off()
  #.tif
  writeRaster(swp.r, paste(filename,"tif",sep="."), format = "GTiff",overwrite=T)
}

###Plot and save standard maps separately
if(standard.maps==T){
  ###Titles and filenames
  title.all = c(paste("a) Early dry season (regular year), monitored SWC = ",round(swc.periods[1,2]),"%",sep=""),
                paste("b) Mid dry season (regular year), monitored SWC = ",round(swc.periods[2,2]),"%",sep=""),
                paste("c) Late dry season (regular year), monitored SWC = ",round(swc.periods[3,2]),"%",sep=""),
                paste("d) Mid dry season (drought year), monitored SWC = ",round(swc.periods[4,2]),"%",sep=""))
  filename.all = c("early_dry_season_regular",
                   "mid_dry_season_regular",
                   "late_dry_season_regular",
                   "mid_dry_season_drought")
  all.maps = list()
  
  ###Create and save all maps            
  for(i in 1:length(swc.periods[,1])){
    print(filename.all[i])
    filename = paste("BCI_SWP_map_",filename.all[i],sep="")
    pred.map.tmp = pred.map
    pred.map.tmp$swc.filled = swc.periods[i,2]
    swp = predict(rf.model, pred.map.tmp)$predicted
    pred.swp = cbind(pred.map.tmp[,c("x","y")],swp)
    pred.swp$swp = ifelse(pred.swp$swp>0,0,pred.swp$swp)
    #.txt
    write.table(pred.swp,paste(filename,"txt",sep="."), row.names=F,sep="\t")
    #.pdf
    swp.r = rasterFromXYZ(pred.swp)
    all.maps[[i]] = swp.r
    pdf(paste(filename,".pdf",sep=""),width = 10, height = 5.69)
    plot(swp.r, col = colors,zlim = zlim,
         smallplot=c(0.88,0.9, 0.3,0.75),
         legend.args=list(text='MPa', side=1, font=1, line=-15.0, cex=1,adj=0))
    elevation.m = matrix(pred.map$elev,ncol = length(x.pred), nrow = length(y.pred),byrow=T)
    maptopo(elevation.m,plotdim=c(1000,500),add=T, interval=2)
    dev.off()
    #.tif
    writeRaster(swp.r, paste(filename,"tif",sep="."), format = "GTiff",overwrite=T)
  }
  
  ###Plot the standards map together
  if(save.figures==T){
    pdf("Figure 3 - All SWP maps.pdf",width = 10,height = 6.12)
    }
  par(mfrow = rep(length(swc.periods[,1])*0.5,2), mar = c(4, 3, 3, 4))
  for(i in 1:length(swc.periods[,1])){
    plot(all.maps[[i]], col = colors,zlim = zlim,legend=F)
    mtext(title.all[i], side=3, line=0.5, at=000,cex=0.80,adj=0)
    if(i==2){plot(all.maps[[i]], col = colors,zlim = zlim, legend.only=T,
                  smallplot=c(0.88,0.9, 0.3,0.75),
                  legend.args=list(text='MPa', side=1, font=1, line=-10.0, cex=1,adj=0))}
    maptopo(elevation.m,plotdim=c(1000,500),add=T, interval=2)
  }
  if(save.figures==T){dev.off()}
}


#.......................................................#
#### Soil sampling locations                         ####
#.......................................................#

#Determine the year and depth of sampling locations for colors
sampling.locations = moisture[,c("date","x","y","depth","location.id")]
sampling.locations$year = as.numeric(format(as.POSIXct(strptime(sampling.locations$date,"%m/%d/%Y",tz="")) ,format = "%Y"))
sampling.x = aggregate(sampling.locations$x,list(sampling.locations$location.id),mean,na.rm=T)
colnames(sampling.x) = c("location.id","x")
sampling.y = aggregate(sampling.locations$y,list(sampling.locations$location.id),mean,na.rm=T)
colnames(sampling.y) = c("location.id","y")
sampling.year = aggregate(sampling.locations$year,list(sampling.locations$location.id),mean,na.rm=T)
colnames(sampling.year) = c("location.id","mean.year")
sampling.depth = aggregate(sampling.locations$depth,list(sampling.locations$location.id),max,na.rm=T)
colnames(sampling.depth) = c("location.id","max.depth")
sampling.aggr = merge(sampling.x,sampling.y)
sampling.aggr = merge(sampling.aggr,sampling.year)
sampling.aggr = merge(sampling.aggr,sampling.depth)
#Plot and save
if(save.figures==T){
  pdf("Figure 1 - Sampling locations.pdf",width = 12, height = 6)
  }
par(mar=c(5.1, 4.1, 4.1, 12.1), xpd=TRUE)
plot(sampling.aggr$x,sampling.aggr$y,col="white",pch=19,xlab=NA,ylab=NA)
sampling.aggr$color = ifelse(sampling.aggr$mean.year=="2015" & sampling.aggr$max.depth<="25","deepskyblue",NA)
sampling.aggr$color = ifelse(sampling.aggr$mean.year=="2015" & sampling.aggr$max.depth>"25","blue",sampling.aggr$color)
sampling.aggr$color = ifelse(sampling.aggr$mean.year=="2016" & sampling.aggr$max.depth<="25","lightgreen",sampling.aggr$color)
sampling.aggr$color = ifelse(sampling.aggr$mean.year=="2016" & sampling.aggr$max.depth>"25","darkgreen",sampling.aggr$color)
sampling.aggr$color = ifelse(sampling.aggr$mean.year<"2016" & sampling.aggr$mean.year>"2015" & sampling.aggr$max.depth<="25","tomato",sampling.aggr$color)
sampling.aggr$color = ifelse(sampling.aggr$mean.year<"2016" & sampling.aggr$mean.year>"2015" & sampling.aggr$max.depth>"25","darkred",sampling.aggr$color)
maptopo(elevation.m,plotdim=c(1000,500),add=T, interval=2,clr="darkgrey")
points(sampling.aggr$x,sampling.aggr$y,bg=sampling.aggr$color,pch=21)
legend("topright", inset=c(-0.25,0), legend=c("2015/16, depth = 15 cm","2015/16, depth > 15 cm",
                                              "2015, depth = 15 cm","2015, depth > 15 cm",
                                              "2016, depth = 15 cm","2016, depth > 15 cm"), pch=rep(21,6),bty="n",
       pt.bg = c("tomato","darkred","deepskyblue","blue","lightgreen","darkgreen"))
if(save.figures==T){dev.off()}


#.......................................................#
#### Soil moisture monitoring and rainfall in 2015-16####
#.......................................................#

#Get monitored SWC in 2015-16
monitored.1516 = monitored.swc.stri
monitored.1516 = monitored.1516[which(monitored.1516$chk_note=="good"),]
monitored.1516 = aggregate(monitored.1516$h2o_by_dry,list(monitored.1516$date),mean,na.rm=T)
colnames(monitored.1516) = c("date","swc")
monitored.1516$date = as.Date(monitored.1516$date, "%d/%m/%Y")
monitored.1516$year = year(monitored.1516$date)
monitored.1516 = monitored.1516[monitored.1516$year %in% c(2015,2016),]
monitored.1516 = monitored.1516[order(monitored.1516$date),]
all.dates = data.frame(date = seq(as.Date(monitored.1516$date[1]), as.Date(monitored.1516$date[length(monitored.1516$date)]), by="days"))
all.dates$date = as.character(all.dates$date)
monitored.1516$date = as.character(monitored.1516$date)
monitored.1516 = merge(all.dates,monitored.1516,all.x=T)
swc.filled = zoo(monitored.1516$swc)
swc.filled = na.approx(swc.filled)
swc.filled = as.numeric(swc.filled)
monitored.1516$swc = swc.filled
monitored.1516$date = as.Date(monitored.1516$date)
#Add rainfall
rain = monitored.rain.stri
rain$date = as.Date(rain$date, "%d/%m/%Y")
rain = rain[,c("date","ra")]
monitored.1516 = merge(monitored.1516,rain,all.x=T)
#Add colors
qu.ra = unique(quantile(monitored.1516$ra,seq(0,1,length.out = 14)))
monitored.1516$class = cut(monitored.1516$ra, breaks = qu.ra, right = FALSE,include.lowest = TRUE)
class.tmp = data.frame(sort(unique(monitored.1516$class)) ,seq(1,8,1))
colnames(class.tmp) = c("class","level")
class.tmp$color = palette(brewer.pal(n = 8, name = "Blues"))
monitored.1516 = merge(monitored.1516,class.tmp,all.x=T)
monitored.1516 = monitored.1516[order(monitored.1516$date),]
monitored.1516$swc.end = c(monitored.1516$swc[2:length(monitored.1516$swc)],NA)
monitored.1516$date.end = c(monitored.1516$date[2:length(monitored.1516$date)],NA)
#Plot
if(save.figures==T){
  pdf("Figure 2 - SWC over time.pdf", width = 10, height = 5)
  par(mar = c(3,6,1,2))
}
ylim = c(20,80)
plot.default(monitored.1516$swc~as.Date(monitored.1516$date),type="l",col = NA,lwd=1,xlab=NA,xaxt="n",  # "black"
     ylab="Soil gravimetric water content\n(% of dry mass)",ylim=ylim,
     panel.first = rect(xleft=as.Date(c("2015/1/1","2015/12/15")), ybottom=rep(ylim[1],2), xlab=NA,
                        xright=as.Date(c("2015/5/1","2016/5/1")), ytop=rep(ylim[2],2), col='lightgrey', border=NA))
segments(x0 = monitored.1516$date, x1 = monitored.1516$date.end, 
         y0 = monitored.1516$swc, y1 = monitored.1516$swc.end,
         col = monitored.1516$color,lwd = 4)
##Add mean of measured SWC in the 50-ha plot per sampling day
measured.swc = moisture[,c("date","swc","period")]
measured.swc$date = as.Date(measured.swc$date, "%m/%d/%Y")
mean.date = aggregate(measured.swc$date,list(measured.swc$period),mean,na.rm=T)
colnames(mean.date) = c("period","meandate")
measured.swc = merge(measured.swc,mean.date)
#Width of boxes should be amount of days per period
width = unique(moisture[,c("date","period")])
width$date = as.Date(width$date, "%m/%d/%Y")
width$day = day(width$date)
min.day = aggregate(width$day, by = list(width$period), FUN = min)
max.day = aggregate(width$day, by = list(width$period), FUN = max)
days = max.day[,2] - min.day[,2] + 1
days = days[c(2,3,1,4)]
#Add boxplots
tmp = boxplot(measured.swc$swc~as.Date(measured.swc$meandate), plot=FALSE)
bxp(tmp, at=sort(unique(as.Date(measured.swc$meandate))),xlim=c(as.Date("2015/1/1"), as.Date("2016/12/31")),outline=F,whisklty = 1,
    boxfill = rgb(.2,.2,.2,alpha=.2),outpch = 19,cex=0.5,axes=F,ylim=ylim,xlab=NA,ylab=NA,add=T,boxwex = days)
axis.Date(side=1, at = seq(as.Date("2015/1/1"), as.Date("2016/12/31"), "months"),labels = FALSE, tcl = -0.25)
axis.Date(side=1, at = seq(as.Date("2015/1/1"), as.Date("2017/1/1"), "years"),labels = TRUE, tcl = -0.75)
#Add legend and save plot
color.legend(as.Date("2016/12/01"),ylim[1],as.Date("2016/12/31"),50,seq(0,180,length.out=8),class.tmp$color,gradient="y",col = rep("white",8),cex=0.001) #
text(as.Date("2016/11/15"),53,"Rainfall (mm)",cex=1)
text(rep(as.Date("2016/10/15"),8),seq(ylim[1]+2,48,length.out=8),class.tmp$class,cex=0.8)
if(save.figures==T){dev.off()}
#Check rainfall on the sampling days.
rain.sampling = rain[which(as.character(rain$date) %in% as.character(monitored.swc.days$date)),]
rain.sampling = cbind(rain.sampling,monitored.swc.days)
rain.sampling[which(rain.sampling$ra>0),]


#.......................................................#
#### Small scale moisture variation                  ####
#.......................................................#

#Calculate SWP (and SWC) distance from centre of the site
centre = small.scale[which(small.scale$distance==0),]
centre = centre[,c("location.id","swp","swc")]
colnames(centre)[2:3] = c("swp.centre","swc.centre")
small.scale = merge(small.scale,centre)
small.scale = small.scale[-which(small.scale$distance==0),]
small.scale$swp.diff = abs(small.scale$swp - small.scale$swp.centre)
small.scale$swc.diff = abs(small.scale$swc - small.scale$swc.centre)
#Create and save plot
theme_new <- theme_set(theme_bw())
theme_new <- theme_update(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_rect(colour = "black"),
                          axis.ticks = element_line(color = "black"),
                          axis.text.y = element_text(size = 12,color = "black"),
                          axis.text.x = element_text(size = 12,color = "black"),
                          axis.title.x = element_text(size = 15,color = "black"),
                          axis.title.y = element_text(size = 15,color = "black"),
                          legend.text = element_blank())
gg.data = small.scale
plot = ggplot(gg.data, aes(factor(gg.data$distance), swp.diff)) +
  geom_violin(fill="grey") +
  xlab("Distance from centre of site (m)") +
  ylab("SWP difference from centre of site (MPa)")
plot
if(save.figures==T){
  ggsave(plot,filename="Figure 4 - Small-scale SWP variation.pdf")
}
theme_set(theme_grey())
##Calculate differences in 10th, 50th and 90th percentiles.
dist.1 = small.scale[which(small.scale$distance==1),]
quantile(dist.1$swp.diff, probs = c(0.10,0.50,0.90),na.rm=T)
dist.2 = small.scale[which(small.scale$distance==2),]
quantile(dist.2$swp.diff, probs = c(0.10,0.50,0.90),na.rm=T)
dist.4 = small.scale[which(small.scale$distance==4),]
quantile(dist.4$swp.diff, probs = c(0.10,0.50,0.90),na.rm=T)


#.......................................................#
#### Sampling size of soils per habitat and soil type####
#.......................................................#

#Determine the percentages of each habitat covered across the plot and by the soil samples.
habitats = bci_habitat  #Available in bciex package
hab.samples = moisture
hab.samples = pick.from.points(data=hab.samples, src = habitats, method = c("nearest.neighbour"), X.name = "x", Y.name = "y", cbind=TRUE)
hab.samples = aggregate(hab.samples$habitat, by = list(hab.samples$habitat), FUN = "length")
colnames(hab.samples) = c("Habitat","n.samples")
hab.samples$perc.samples = hab.samples$n.samples / sum(hab.samples$n.samples) * 100
hab.area = aggregate(habitats$habitat, by = list(habitats$habitat), FUN = "length")
colnames(hab.area) = c("Habitat","perc.area")
hab.area$perc.area = hab.area$perc.area / sum(hab.area$perc.area) * 100
hab.all = merge(hab.samples,hab.area)
hab.all$Habitat = c("High plateau","Low plateau","Mixed","Slope","Stream","Swamp","Young")
hab.all$Habitat.soil = "Habitat"
colnames(hab.all)[which(colnames(hab.all)=="Habitat")] <- "Category"
hab.all = hab.all[,c("Habitat.soil","Category","n.samples","perc.samples","perc.area")]
##Do the same for soil types
soil.samples = moisture
soil.samples = pick.from.points(data=soil.samples, src = soiltype, method = c("nearest.neighbour"), X.name = "x", Y.name = "y", cbind=TRUE)
soil.samples = aggregate(soil.samples$soiltype, by = list(soil.samples$soiltype), FUN = "length")
colnames(soil.samples) = c("Soil_type","n.samples")
soil.samples$perc.samples = soil.samples$n.samples / sum(soil.samples$n.samples) * 100
soil.area = aggregate(soiltype$soiltype, by = list(soiltype$soiltype), FUN = "length")
colnames(soil.area) = c("Soil_type","perc.area")
soil.area$perc.area = soil.area$perc.area / sum(soil.area$perc.area) * 100
soil.all = merge(soil.samples,soil.area)
soil.all$Soil_type = c("AVA","Fairchild","Marron","Swamp")
soil.all$Habitat.soil = "Soil type"
colnames(soil.all)[which(colnames(soil.all)=="Soil_type")] <- "Category"
soil.all = soil.all[,c("Habitat.soil","Category","n.samples","perc.samples","perc.area")]
#Merge and save table
soil.hab.all = rbind(hab.all,soil.all)
filename = "Table 2 - Samples and area.xlsx"
if(save.figures==T){write.xlsx(soil.hab.all, filename,row.names=F)}

