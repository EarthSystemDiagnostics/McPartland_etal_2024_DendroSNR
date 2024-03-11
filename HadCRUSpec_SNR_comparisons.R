library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(extrafont)
library(rgdal)
library(lattice)
library(maps)
library(fields)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(maptools)
library(reshape2)
library(stars)
library(gridExtra)
library(raster)
library(ggridges)
library(ggpubr)
library(cowplot)
library(zoo)
library(tmap)
library(ncdf4.helpers)
library(ncdf4)
require(lemon)
require(ggforce)
#install.packages("PaleoSpec")
require(PaleoSpec)
source("./FilterSpec.R")
########interpolation function 

SpecApprox<-function(spec, xout=NULL,...){
  if(is.null(xout)){
    #frq.bnds<-log10(range(spec$freq))
    frq.bnds<-c(1/1000,1/2) #bin output by regular frequency bands
    xout<-seq(frq.bnds[1],frq.bnds[2], 0.001)
  }
  int.spec<-approx(x = spec$freq,y = spec$spec,xout = xout,...)
  xnt.spec<-approx(x=spec$freq, y=spec$dof, xout = xout,...)
  #znt.spec<-approx(x=spec$freq, y=spec$Scales, xout = xout,...)
  rtrn.spec<-list(freq=int.spec$x,
                  spec=int.spec$y,
                  dofs=xnt.spec$y)
  #scales=znt.spec$y)
  class(rtrn.spec)<-"spec"
  return(rtrn.spec)
}

####load datasets
#source global temps dataset
global_temps <- nc_open("HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean_INFILLED.nc", write=FALSE, readunlim=TRUE, verbose=FALSE, 
                        auto_GMT=TRUE, suppress_dimvals=FALSE)
points<-read.csv("meta_sites_inSNR.csv")
snrs<-read.csv("bootstrapped_curves.csv")


##########calculate mean summer temp and extract grid cells from hadcru netcdf file
lon <- ncvar_get(global_temps,"longitude")
lat <- ncvar_get(global_temps,"latitude")

# (Starts 01/01/1850)
time <- ncvar_get(global_temps,"time")

#use this to identify the first timestep
time_d <- as.Date(time, format="%j", origin=as.Date("1850-01-01"))
time_c<-as.character(time_d)

#use ncvar_get to get the temps data
temps_array <- ncvar_get(global_temps, "tas_mean")

#find locations for the starting lat lon and time (only if sub-setting by AOI or time slice, here we take full global record) 
start_n <- which(time_c == "1850-01-16")
start<-which(time == 15.5)
start_lon <- which(lon == -177.5)
start_lat <- which(lat == -87.5)

start<-c(1,1,1) #first two numbers refer to the lon and lat, the third
#refers to the year that you want to cell to start at (i.e. the first month of the first year)

#first two numbers refer to the lengths of lat and lan, and the third refers to 
#the number of years you want to incorporate in the analysis
#cut it off in dec of 2021
count=c(72, 36, 2064)


#subset data by bounds
temps_subset <- ncvar_get(global_temps, "tas_mean", 
                          start=start,
                          count=count)

#create year and month vector to sort and filter and average
Year <- rep(seq(1850,2021),each=12)
Month <- rep(seq(1,12),times=172)

#define dimensions and transpose temps data
dim(temps_subset)          
dim(temps_subset) <- c(72*36,2064)
temps_subset <- t(temps_subset)

#bind temps and temporal data
temps_merge <- cbind(Year,Month, temps_subset)

#convert to data frame
temps_merge <- as.data.frame(temps_merge)

#filter summer months
temps_filter <- temps_merge %>%
  filter(Month == 6 | Month == 7 | Month == 8)

#average temps by year
temps_mean <- temps_filter %>%
  group_by(Year) %>%
  summarise_all(list(mean), na.rm = T)


#remove month because it doesn't matter now
temps_mean$Month<-NULL
temps_mean[temps_mean == "NaN"] <-NA

#convert to list to iterate without dataframe constraints
temps_list<-as.list(temps_mean)

#save and remove year from list
years<-temps_list$Year
temps_list$Year<-NULL


##north america bounding box (-178.2, 6.6, -49, 83.3)
#turn into maxtrix and convert to raster with appropriate dimensions. 
temps_mean_matrix <- matrix(unlist(temps_list), nrow = 2592, byrow = T)
temps_brick <- brick(nrow=36, ncol=72, nl=172)
values(temps_brick)<-temps_mean_matrix
extent(temps_brick) <- c(-177.5, 177.5, -87.5 , 87.5)
temps_brick<-flip(temps_brick)
#visualize random year in brick to make sure everything is lined up correctly
image(temps_brick$layer.170)

#convert lat lon from metadata to spatial object 
lonlat<-points[,5:4]
plot(lonlat)
pts<-SpatialPoints(lonlat, proj4string = CRS("+proj=longlat +datum=WGS84"))
plot(pts)
lonlat$index_cells<-rownames(lat)

#extract points from raster. 
cru_extraction<-extract(temps_brick, pts, cellnumbers = T)


#plot original HADCRU points 
original_cru_coords<-coordinates(temps_brick)[cru_extraction[,1],]
plot(original_cru_coords)

#make the extracted time series into a dataframe and then a list
temps_extract<-as.data.frame(cru_extraction)

#index cells save the spatial references from the raster
index_cells<-temps_extract$cells

#remove the index cells and add colnames as years
temps_extract$cells<-NULL
colnames(temps_extract)<-years


#site level and HadCRU coordinates
original_cru_coords<-as.data.frame(original_cru_coords)
original_cru_coords$pages_lat<-lonlat$lat
original_cru_coords$pages_lon<-lonlat$lon
original_cru_coords$name<-paste0(tolower(points$code))

#take temps extract, make into list and add names from sites
temps_extract<-as.data.frame(t(temps_extract))
temps_extract_list<-as.list(temps_extract)
names(temps_extract_list)<-paste0(points$code)
temps_extract_list$Year<-years
hadcru<-temps_extract_list

meta<-points

##clear environment for 
rm(list=setdiff(ls(), c("meta","hadcru","SpecApprox","snrs","FilterSpec","FilterSpecLog")))

###################process HADCRU data to interpolate
hadcru_years<-hadcru$Year
hadcru$Year<-NULL

#make sure the correct sites are represented
diff<-setdiff(names(hadcru), meta$code)
hadcru<-hadcru[!names(hadcru) %in% diff]

#filter raw values from hadcru raw data to remove any non-contiguous regions
#this looses it's association with year, but then the series are converted to spectra so it's not really an issue

hadcru_filter<-list()

for(i in hadcru) {
  filt<-if(all(is.na(i))) NA else((na.contiguous(i)))
  hadcru_filter[[length(hadcru_filter)+1]]<-unlist(filt)
}

names<-names(hadcru)
hadcru_names<-paste0(names, "-hadcru")

names(hadcru_filter)<-hadcru_names

#bind all dataframes together
clim<-c(hadcru_filter)

clim<-clim[!duplicated(clim)]

#perform the spectral analysis of model and instrumental data
clim_specs<-list()

for(i in 1:length(clim)) {
  tmp<-clim[[i]]
  name <- names(clim[i])
  print(name)
  tmp_ts<-ts(tmp, deltat = 1)
  spectrum <- SpecMTM(tmp_ts)
  spc.apprx<-SpecApprox(spectrum)
  name_split<-strsplit(name, "-")
  clim_specs[[paste0(name)]] <- spc.apprx
}  

###make data frame of climate data
clim_df<-setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("freq","spec","dofs", "data_type","proxy","name"))

for(i in 1:length(clim_specs)) {
  freq <- clim_specs[[i]]$freq
  spec <- clim_specs[[i]]$spec
  dofs<-clim_specs[[i]]$dofs
  name <- names(clim_specs[i])
  sitename<-strsplit(name, '-')[[1]][1]
  data_type<-strsplit(name, '-')[[1]][2]
  spec_bind<-cbind(freq, spec, dofs, data_type, sitename)
  spec_bind<-as.data.frame(spec_bind)
  clim_df<-rbind(clim_df, spec_bind)
}

clim_df$freq<-as.numeric(as.character(clim_df$freq))
clim_df$spec<-as.numeric(as.character(clim_df$spec))
clim_df$dofs<-as.numeric(as.character(clim_df$dofs))

#visualze all specs
test_plot<-ggplot(clim_df, aes(x=freq, y = spec, group = sitename, color = data_type))+
  geom_line()+
  scale_x_continuous(trans = c('log10'))+
  scale_y_continuous(trans = c('log10'))+
  facet_grid(.~data_type)
test_plot

############### 
#take mean and add DOF of 1
clim_sums <- clim_df %>%
  group_by(data_type, freq) %>%
  summarise(spec = mean(spec, na.rm = T),
            dofs = 1)%>%
  filter(spec != "NaN")


theme1<-theme(#legend.key = element_blank(),
  text = element_text(size = 24),
  title = element_text(size = 24),
  strip.text.x = element_text(size = 20),
  strip.background.x = element_blank(),
  strip.text.y = element_text(size = 20))

######filter less than 100 yr timescales

hadcru<-clim_sums%>%
  filter(data_type == "hadcru" & freq >= 0.01)%>%
  filter(spec != "NaN")

#smooth spectra
spec_list<-list(hadcru)

for(i in 1:length(spec_list)){
  df<-spec_list[[i]]
  ser<-list(df$spec, df$freq, df$dofs)
  names(ser)<-c("spec","freq","dof")
  class(ser)<-"spectrum"
  smooth<-FilterSpec(ser, spans = c(3,5))
  smooth<-FilterSpecLog(smooth, df = 0.5)
  df$spec<-smooth$spec
  df$upper<-smooth$lim.1
  df$lower<-smooth$lim.2
  spec_list[[i]]<-df
}

specs_smoothed<-do.call(rbind, spec_list)

#filter out low-frequencies from SNR data
snrs<-snrs%>%
  filter(proxy == "MXD" | freq >= 0.01)%>%
  filter(proxy == "TRW" | freq >= 0.01)

snrs<-snrs%>%
  filter(!is.na(mean_signal_curve))

##additional smoothing on mean signal and raw spectra for MXD and TRW
mxd<-snrs%>%
  filter(proxy == "MXD")
mxd<-list(mxd$freq, mxd$mean_signal_curve)
names(mxd)<-c("freq","spec")
mxd$dofs<-rep(5, length(mxd$freq))
class(mxd)<-"spectrum"
mxd<-FilterSpec(mxd, spans = c(3,5))
mxd<-FilterSpecLog(mxd, df = 0.5)
snrs$mean_signal_curve[snrs$proxy == "MXD"]<-mxd$spec

mxd<-snrs%>%
  filter(proxy == "MXD")
mxd<-list(mxd$freq, mxd$mean_raw_curve)
names(mxd)<-c("freq","spec")
mxd$dofs<-rep(5, length(mxd$freq))
class(mxd)<-"spectrum"
mxd<-FilterSpec(mxd, spans = c(3,5))
mxd<-FilterSpecLog(mxd, df = 0.5)
snrs$mean_raw_curve[snrs$proxy == "MXD"]<-mxd$spec

trw<-snrs%>%
  filter(proxy == "TRW")
trw<-list(trw$freq, trw$mean_signal_curve)
names(trw)<-c("freq","spec")
trw$dofs<-rep(5, length(trw$freq))
class(trw)<-"spectrum"
trw<-FilterSpec(trw, spans = c(3,5))
trw<-FilterSpecLog(trw, df = 0.5)
trw$mean_signal_curve[snrs$proxy == "TRW"]<-trw$spec

trw<-snrs%>%
  filter(proxy == "TRW")
trw<-list(trw$freq, trw$mean_raw_curve)
names(trw)<-c("freq","spec")
trw$dofs<-rep(5, length(trw$freq))
class(trw)<-"spectrum"
trw<-FilterSpec(trw, spans = c(3,5))
trw<-FilterSpecLog(trw, df = 0.5)
trw$mean_raw_curve[snrs$proxy == "TRW"]<-trw$spec


#add confidence intervals from bootstrap to real estimates through multiplication 
snrs$signal_0.9<-snrs$mean_signal_curve*snrs$signal_0.9
snrs$signal_0.1<-snrs$mean_signal_curve*snrs$signal_0.1

#selected series for rescaling and plotting 
snrs<-snrs %>%
  dplyr::select(freq, mean_signal_curve, signal_0.1, signal_0.9, proxy)

#combine hadcrut data
specs_smoothed<-cbind.data.frame(specs_smoothed$freq, specs_smoothed$spec, specs_smoothed$lower, specs_smoothed$upper, specs_smoothed$data_type)

#filter out raw series and interannual timescales
snrs_raw<-snrs
snrs_raw<-snrs_raw%>%
  filter(proxy == "TRW" | freq >= 0.0125)

snrs_raw<-snrs_raw%>%
  dplyr::select(freq, mean_raw_curve, rawLowerCI, rawUpperCI, proxy)

snrs_raw$proxy[snrs_raw$proxy == "MXD"]<-"MXD_raw"
snrs_raw$proxy[snrs_raw$proxy == "TRW"]<-"TRW_raw"

colnames(snrs_raw)<-c("freq","spec","lower","upper","data_type")
colnames(snrs)<-c("freq","spec","lower","upper","data_type")
colnames(specs_smoothed)<-c("freq","spec","lower","upper","data_type")

#combine all datatypes 
all_specs<-rbind(snrs, specs_smoothed, snrs_raw)
all_specs$lower[all_specs$data_type == "TRW_raw" | all_specs$data_type == "MXD_raw"]<-NA
all_specs$upper[all_specs$data_type == "TRW_raw" | all_specs$data_type == "MXD_raw"]<-NA

#mean across all spectral estimates (0-8 years)
mean_all<-mean(all_specs$spec[all_specs$freq <= 0.125 & all_specs$freq >= 0.01])

spec_high<-all_specs %>%
  group_by(data_type)%>%
  filter(freq <= 0.125 & freq >= 0.01) %>%
  summarise(mean_high = mean(spec),
            scaled = 1/(mean_all/mean_high))

#normalise all indiviudal estimates by the global mean
all_specs$spec[all_specs$data_type == "TRW"]<-all_specs$spec[all_specs$data_type == "TRW"]/spec_high$scaled[spec_high$data_type == "TRW"]
all_specs$spec[all_specs$data_type == "MXD"]<-all_specs$spec[all_specs$data_type == "MXD"]/spec_high$scaled[spec_high$data_type == "MXD"]
all_specs$spec[all_specs$data_type == "hadcru"]<-all_specs$spec[all_specs$data_type == "hadcru"]/spec_high$scaled[spec_high$data_type == "hadcru"]


all_specs$upper[all_specs$data_type == "TRW"]<-all_specs$upper[all_specs$data_type == "TRW"]/spec_high$scaled[spec_high$data_type == "TRW"]
all_specs$upper[all_specs$data_type == "MXD"]<-all_specs$upper[all_specs$data_type == "MXD"]/spec_high$scaled[spec_high$data_type == "MXD"]
all_specs$upper[all_specs$data_type == "hadcru"]<-all_specs$upper[all_specs$data_type == "hadcru"]/spec_high$scaled[spec_high$data_type == "hadcru"]


all_specs$lower[all_specs$data_type == "TRW"]<-all_specs$lower[all_specs$data_type == "TRW"]/spec_high$scaled[spec_high$data_type == "TRW"]
all_specs$lower[all_specs$data_type == "MXD"]<-all_specs$lower[all_specs$data_type == "MXD"]/spec_high$scaled[spec_high$data_type == "MXD"]
all_specs$lower[all_specs$data_type == "hadcru"]<-all_specs$lower[all_specs$data_type == "hadcru"]/spec_high$scaled[spec_high$data_type == "hadcru"]

all_specs$spec[all_specs$data_type == "TRW_raw"]<-all_specs$spec[all_specs$data_type == "TRW_raw"]/spec_high$scaled[spec_high$data_type == "TRW_raw"]
all_specs$spec[all_specs$data_type == "MXD_raw"]<-all_specs$spec[all_specs$data_type == "MXD_raw"]/spec_high$scaled[spec_high$data_type == "MXD_raw"]


p2<-ggplot()+
  theme_bw()+
  geom_line(data = all_specs, aes(x=freq, y=spec, group = data_type, color = data_type, linetype = data_type), size = 1.5)+
  #geom_line(data = snrs_raw, aes(x=freq, y = mean_raw_curve, color = proxy))+
  geom_ribbon(data = all_specs, aes(x = freq, ymin = lower , ymax = upper, group = data_type, fill = data_type), alpha = 0.5, size = 1)+
  #geom_ribbon(data = snrs_raw, aes(x = freq, ymin = mean_raw_curve - rawLowerCI , ymax = mean_raw_curve+rawUpperCI, group = proxy, fill = proxy), alpha = 0.5, size = 1)+
  theme(panel.border=element_blank(), axis.line=element_line())+
  ylab("PSD")+
  xlab("Timescale")+
  scale_fill_manual(name = "Data Type", breaks = c("TRW","TRW_raw","MXD","MXD_raw","hadcru"),
                    labels = c("Tree-ring width\n(corrected)","Tree-ring width (raw)","Tree-ring density\n(corrected)","Tree-ring density (raw)", "HadCRUT Summer\nTemp Anomalies"),
                    values = c("#008837",NA, "#7b3294",NA, "darkgrey"), na.value = NA)+
  scale_color_manual(name = "Data Type", breaks = c("TRW","TRW_raw","MXD","MXD_raw","hadcru"),
                     labels = c("Tree-ring width\n(corrected)","Tree-ring width (raw)","Tree-ring density\n(corrected)","Tree-ring density (raw)", "HadCRUT Summer\nTemp Anomalies"),
                     values = c("#008837","#a6dba0", "#7b3294","#c2a5cf", "darkgrey"), na.value = NA)+
  scale_linetype_manual(name = "Data Type", breaks = c("TRW","TRW_raw","MXD","MXD_raw","hadcru"),
                        labels = c("Tree-ring width\n(corrected)","Tree-ring width (raw)","Tree-ring density\n(corrected)","Tree-ring density (raw)", "HadCRUT Summer\nTemp Anomalies"),
                        values = c("solid","dashed", "solid","dashed", "solid"))+
  scale_x_continuous(trans = c('log10'), limits = c(0.01,0.50),
                     breaks = c(0.01,0.033, 0.1, 0.50),
                     labels = c('100',"30","10","2"), 
                     sec.axis = sec_axis((~.), name = "Frequency (1/years)", 
                                         breaks = c(0.01,0.033, 0.1, 0.5), 
                                         labels = c("0.01","0.03","0.1","0.5")))+
  scale_y_continuous(trans=c('log10'), limits = c(0.03,2.8))+
  coord_capped_cart(bottom="both", left="both", top = "both")+
  ggtitle("(a)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        legend.position = c(0.01,0.17),
        legend.justification = c(0),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.direction="horizontal",
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0),
        strip.text.x = element_text(size = 12),
        strip.background.x = element_blank(),
        strip.text.y = element_text(size = 12))+
  guides(fill=guide_legend(nrow=3,byrow=TRUE), label.position = "right")
p2

