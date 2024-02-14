#load packages
require(dplyr)
require(tidyr)
require(ggplot2)
require(extrafont) 
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(scales)
require(graphics)
require(utils)
require(lemon)
require(ggforce)
#devtools::install_github("EarthSystemDiagnostics/paleospec")
require(PaleoSpec)
require(sp)
require(rgdal)
library(geosphere)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(maptools)
library(reshape2)
library(stars)
library(raster)
library(ggridges)
library(ggpubr)
library(zoo)
library(tmap)
library(ncdf4.helpers)
library(ncdf4)
source("./FilterSpec.R")
####load approximation function
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

##load pages database from NOAA website
tree<-load("/Users/mamcpa001/Downloads/Pages2kTemperature2_1_2.RData")

#filter tree-ring records
variableName <- as.character(sapply(D, "[[", "archiveType"))
tree<-which(variableName == "tree")
D <- D[tree]

rm(sTS)
rm(TS)


#append lat and lon info to each chronology, and to each site
for(i in 1:length(D)){
  obj<-D[[i]]
  df_name<-names(D)[[i]]
  name<-names(obj)
  lat<-obj$geo$latitude
  lon<-obj$geo$longitude
  obj$paleoData<-lapply(obj$paleoData, function(x) {x$lon <- lon; x})
  obj$paleoData<-lapply(obj$paleoData, function(x) {x$lat <- lat; x})
  obj$lat<-lat
  obj$lon<-lon
  D[[i]]<-obj
}

#filter sites in northern hemisphere
lat_var <- as.numeric(sapply(D, "[[", "lat"))
lat<-which(lat_var > 0)
D <- D[lat]

#optional filter sites in western hemisphere  for only N. American sites
#lon_var <- as.numeric(sapply(D, "[[", "lon"))
#lon<-which(lon_var > 0)
#D <- D[lon]

#simplfy database to only paleo information and unlist to access individual chronologies
paleoData<-sapply(D, "[[", "paleoData")
paleoData<-unlist(paleoData, recursive = F)
rm(D)

#filter out specific proxy types
#measurementTable$MXD is unitless (standardized)
#deltaDensity (only a few crns) is unitless (standardized) within measurementTable$density

density<-paleoData[-which(sapply(paleoData, function(x) (is.null(x$measurementTable[[1]]$density))))]
mxd<-paleoData[-which(sapply(paleoData, function(x) (is.null(x$measurementTable[[1]]$MXD))))]

density_names<-names(density)
mxd_names<-names(mxd)

diff<-setdiff(density_names, mxd_names)

#keep only the standardized delta density sites
density<-density[1:5]

#select tree ring sites
trw<-paleoData[-which(sapply(paleoData, function(x) (is.null(x$measurementTable[[1]]$trsgi))))]



# for each proxy type, extract each chronology and create new object that contains 
# standardized values, year, sitename, lat, lon, proxy, and (if available) sample depth, rbar, and eps
mxd_list<-list()

for(i in 1:length(mxd)){
  obj<-mxd[[i]]
  sitename<-obj$measurementTable[[1]]$MXD$TSid
  proxy<-obj$measurementTable[[1]]$MXD$proxy
  lat<-obj$lat
  lon<-obj$lon
  values<-obj$measurementTable[[1]]$MXD$values
  samplecount<-obj$measurementTable[[1]]$sampleCount$values
  rbar<-obj$measurementTable[[1]]$RBar$values
  eps<-obj$measurementTable[[1]]$EPS$values
  year<-obj$measurementTable[[1]]$year$values
  new<-list(sitename, proxy, lat, lon, values, samplecount,rbar,eps,year)
  names(new)<-c("code","proxy","lat","lon","rwi","samplecount","rbar","eps","year")
  new<-list(new)
  names(new)<-sitename
  mxd_list<-c(mxd_list, new)
}

trw_list<-list()

for(i in 1:length(trw)){
  obj<-trw[[i]]
  sitename<-obj$measurementTable[[1]]$trsgi$TSid
  proxy<-obj$measurementTable[[1]]$trsgi$proxy
  lat<-obj$lat
  lon<-obj$lon
  values<-obj$measurementTable[[1]]$trsgi$values
  samplecount<-obj$measurementTable[[1]]$sampleCount$values
  rbar<-obj$measurementTable[[1]]$RBar$values
  eps<-obj$measurementTable[[1]]$EPS$values
  year<-obj$measurementTable[[1]]$year$values
  new<-list(sitename, proxy, lat, lon, values, samplecount,rbar,eps,year)
  names(new)<-c("code","proxy","lat","lon","rwi","samplecount","rbar","eps","year")
  new<-list(new)
  names(new)<-sitename
  trw_list<-c(trw_list, new)
}


density_list<-list()

for(i in 1:length(density)){
  obj<-density[[i]]
  sitename<-obj$measurementTable[[1]]$density$TSid
  proxy<-obj$measurementTable[[1]]$density$proxy
  lat<-obj$lat
  lon<-obj$lon
  values<-obj$measurementTable[[1]]$density$values
  samplecount<-obj$measurementTable[[1]]$sampleCount$values
  rbar<-obj$measurementTable[[1]]$RBar$values
  eps<-obj$measurementTable[[1]]$EPS$values
  year<-obj$measurementTable[[1]]$year$values
  new<-list(sitename, proxy, lat, lon, values, samplecount,rbar,eps,year)
  names(new)<-c("code","proxy","lat","lon","rwi","samplecount","rbar","eps","year")
  new<-list(new)
  names(new)<-sitename
  density_list<-c(density_list, new)
}


#combine all in one list
pages<-c(mxd_list, density_list, trw_list)


#create metadata file for new database
lat<-sapply(pages, "[[", "lat")
lon<-sapply(pages, "[[", "lon")
name<-names(pages)
proxy<-sapply(pages, "[[", "proxy")
year<-sapply(pages, "[[", "year")
min_year<-sapply(year, min)
max_year<-sapply(year, max)


meta_all<-cbind.data.frame(name, proxy, lat, lon, min_year, max_year)
colnames(meta_all)<-c("code","proxy","lat","lon", "startyear","endyear")

rm(paleoData)
rm(density_list)
rm(mxd_list)
rm(trw_list)


#write.csv(meta, file = "meta_pages17.csv")
#saveRDS(pages, file = "pages17_extract.rds")

#############################
#define spatial clusters of sites for TRW and density data separately 
##THIS SECTION SHOULD BE RUN TWICE, ONCE FOR TRW AND ONCE FOR DENSITY SITES

#run for tree ring sites only
#meta<-meta_all%>%
#  filter(proxy == "TRW")

#run for density sites
meta<-meta_all%>%
  filter(proxy=="MXD" | proxy=="delta Density")

location<-as.data.frame(as.character(meta$code))
location$latitude<-as.numeric(meta$lat)
location$longitude<-as.numeric(meta$lon)
colnames(location)[1]<-"code"
id<-names(pages)

#calculate the distances between all sites and return large matrix
m <- distm(location[3:2], location[3:2], fun = distHaversine)

#remove half of observations which are redundant
diag(m) <- NA
m <- m * lower.tri(m)
m[m == 0] <- NA

#turn into dataframe with testament on the cols and rows
names<-c(paste(location$code))
m<-as.data.frame(m)
colnames(m)<-names
rownames(m)<-names

#meters to kilometers
m<-m/1000

#all points within or equal to 250km of each other (true/false binary)
m<-m<=250
m<-as.data.frame(m)

m$names<-rownames(m)

#gather up the true/false matches
m<-m %>%
  gather("names2" , "match", 1:length(m)-1)

#remove empty cells which are their own match, results in all sites that are within 250 km of each other
m<-m %>%
  filter(!is.na(match))

#remove non-matches
m<-m %>%
  filter(match != FALSE)

m$match<-NULL

#create a list that contains each unique site with all it's matches as a list
sites_list <-m %>%
  pivot_wider(names_from = names, values_from = names2, values_fn = list)

#make into a list
sites_list<-as.list(sites_list)
#un-nest one level
sites_list <- unlist(sites_list,recursive=FALSE)

#remove any clusters that contain less than or equal to 3 sites (including the pivot site, this eliminates all clusters of less than 4)
sites_list<-sites_list[sapply(sites_list, function(x) all(length(x)>=2))]

names<-names(sites_list)

#add the pivot site to the list, and re-name each cluster
for(i in 1:length(sites_list)){
  sites_list[[i]]<-c(sites_list[[i]], names[i])
}

#remove pivot site names
names(sites_list)<-NULL

#assign each cluster a number
clusters<-paste0(seq_along(1:length(sites_list)))
clusters<-as.numeric(clusters)

names(sites_list)<-clusters

clusters<-sites_list

#go into into the processed pages database and extract all sitenames within each cluster and contain them in a listed object
all_clusters<-list()

for(i in 1:length(clusters)){
  clust<-clusters[i]
  cluster_name<-names(clust)
  clust<-unlist(clust, recursive = FALSE)
  diff<-setdiff(id, clust)
  new_clust<-pages[!names(pages) %in% diff]
  new_clust<-lapply(new_clust, function(x)  {x$id <- paste(x$name, x$proxy);x})
  new_clust<-list(new_clust)
  names(new_clust)<-cluster_name
  all_clusters<-c(all_clusters, new_clust)
}


mxd<-all_clusters
#trw<-all_clusters

#####################combine trw and mxd into one data object for the main analysis
all_clusters<-c(mxd, trw)
########save names of sites included in clusters to extract and analysis along with instrumental record.
# sites_in_clusters<-unlist(all_clusters, recursive=F)
# sites_in_clusters<-lapply(sites_in_clusters, function(x) x$code)
# sites_in_clusters<-unlist(sites_in_clusters)
# sites_in_clusters<-unique(sites_in_clusters)
# meta_in_SNR<-meta_all%>%
#   filter(code %in% sites_in_clusters)
# write.csv(meta_in_SNR, file = "meta_in_SNR.csv")

numbers<-rep(1:length(all_clusters))
proxy<-lapply(all_clusters, function(x) x[[1]]$proxy[1])
all_names<-lapply(all_clusters, function(x) names(x))

names(all_clusters)<-paste0(numbers, "_",proxy)


#transform each series into a single spectral estimate, first making sure each series is the same length
clusters_specs<-list()

for(i in 1:length(all_clusters)) {
  site<-all_clusters[i]
  cluster_name<-names(site)
  site<-unlist(site,recursive=FALSE)
  spec_list<-list()
  min_year<-lapply(site, function(x) min(x$year))
  min_year<-max(unlist(min_year))
  max_year<-lapply(site, function(x) max(x$year))
  max_year<-min(unlist(max_year))
  for(j in 1:length(site)) {
    rwi <- site[[j]]$rwi
    name <- site[[j]]$code
    proxy <- site[[j]]$proxy[1]
    year<-site[[j]]$year
    rwi<-as.data.frame(cbind(rwi, year))
    rwi<-rwi%>%
      filter(year >= min_year)%>%
      filter(year <= max_year)
    rwi<-ts(rwi$rwi, deltat = 1)
    spectrum<-SpecMTM(rwi)
    spectrum$freq<-spectrum$freq[-(1:2)]
    spectrum$spec<-spectrum$spec[-(1:2)]
    spectrum$dof<-spectrum$dof[-(1:2)]
    spec_list[[paste0(name, "_" , proxy)]] <- spectrum
    spec_list[[j]]$name<-name
  }
  spec_list<-list(spec_list)
  print(cluster_name)
  names(spec_list)<-cluster_name
  clusters_specs<-c(clusters_specs, spec_list)
}

####for each cluster, calculate the mean of all spectra
##list with each individual spectral estimate
df_list<-list()

#list with mean curves, along with the sum of the degrees of freedom and number of spectral estimates in the mean curve calculation
spec_long_list<-list()

for(j in 1:length(clusters_specs)){
  spec_list<-clusters_specs[j]
  cluster_name<-names(spec_list)
  spec_list<-unlist(spec_list,recursive=FALSE)
  spec_df<-setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("freq","spec" ,"dofs","name"))
  for(i in 1:length(spec_list)) {
    freq <- spec_list[[i]]$freq
    spec <- spec_list[[i]]$spec
    dof <- spec_list[[i]]$dof
    n <- length(spec_list)
    name <- replicate(length(freq), names(spec_list[i]))
    spec_bind<-cbind(freq, spec, dof, name)
    spec_bind<-as.data.frame(spec_bind)
    spec_df<-rbind(spec_df, spec_bind)
    spec_df$freq<-as.numeric(spec_df$freq)
    spec_df$spec<-as.numeric(spec_df$spec)
    spec_df$dof<-as.numeric(spec_df$dof)
    spec_df <- spec_df%>%
      filter(!is.na(spec))
    spec_sum <- spec_df %>%
      group_by(freq) %>%
      summarise(mean_of_all_spectra = mean(spec, na.rm = T),
                sum_dofs = mean(dof, na.rm = T),
                n.spec = n())
    spec_sum$n.spec<-n
    spec_sum<-as.data.frame(spec_sum)
    spec_sum<-spec_sum[,c("freq","mean_of_all_spectra", "sum_dofs", "n.spec")]
    spec_smooth<-list(spec_sum$freq, spec_sum$mean_of_all_spectra, spec_sum$sum_dofs)
    names(spec_smooth)<-c("freq","spec","dof")
    class(spec_smooth)<-"spectrum"
    spec_smooth<-FilterSpec(spec_smooth, spans = c(3,5))
    spec_smooth<-FilterSpecLog(spec_smooth, df = 0.01)
    spec_smooth<-AddConfInterval(spec_smooth)
    spec_sum$mean_of_all_spectra<-spec_smooth$spec
    spec_sum$upperCI<-spec_smooth$lim.1
    spec_sum$lowerCI<-spec_smooth$lim.2
    spec_sum$name<-cluster_name
    spec_long<-spec_sum %>%
      gather("type", "val", 2:4)
    spec_long<-spec_long[,c("freq","val","type")]
    final<-list(spec_df)
    spec_long<-list(spec_long)
    names(final)<-cluster_name
    names(spec_long)<-cluster_name
  }
  df_list<-c(df_list, final)
  spec_long_list<-c(spec_long_list, spec_long)
}


#format all sites in cluster and truncate to overlapping period
format_list<-list()
#create a list with min/max years for each cluster, along with other relevant data (if available) like sample density
#also calculate the correlation between all series in the cluster
minmaxyears<-list()

for(j in 1:length(all_clusters)){
  sites<-all_clusters[j]
  cluster_name<-names(sites)
  sites<-unlist(sites,recursive=FALSE)
  min_year<-lapply(sites, function(x) min(x$year))
  min_year<-max(unlist(min_year))
  max_year<-lapply(sites, function(x) max(x$year))
  max_year<-min(unlist(max_year))
  mean_ncor<-lapply(sites, function(x) mean(x$samplecount))
  mean_ncor<-mean(unlist(mean_ncor), na.rm=T)
  mean_rbar<-lapply(sites, function(x) mean(x$rbar))
  mean_rbar<-mean(unlist(mean_rbar), na.rm = T)
  mean_eps<-lapply(sites, function(x) mean(x$eps))
  mean_eps<-mean(unlist(mean_eps), na.rm = T)
  dfs<-lapply(sites, function(x) cbind(x$rwi, replicate(length(x$year), proxy), x$year, x$code))
  dfs<-lapply(dfs, function(x) as.data.frame(x))
  colnames<-c("rwi","proxy", "year", "code")
  dfs<-lapply(dfs, setNames, nm = colnames)
  dfs <- lapply(dfs, function(x) {x$rwi <- as.numeric(x$rwi);x})
  dfs <- lapply(dfs, function(x) {x$year <- as.numeric(x$year);x})
  dfs<-lapply(dfs, function(x) filter(x, year <= max_year))
  dfs<-lapply(dfs, function(x) filter(x, year >= min_year))
  cor_df<-lapply(dfs, function(x) cbind.data.frame(x$rwi))
  cor_df<-do.call(cbind, cor_df)
  cor_df<-round(cor(cor_df[, 1:ncol(cor_df)], use="pair"),2)
  diag(cor_df) <- NA
  cor_df <- cor_df*lower.tri(cor_df)
  cor_df[cor_df == 0] <- NA
  cor<-mean(cor_df, na.rm = T)
  minmax<-list(c(min_year, max_year, cor, mean_ncor, mean_rbar, mean_eps))
  dfs<-list(dfs)
  names(dfs)<-cluster_name 
  format_list<-c(format_list, dfs)
  minmaxyears<-c(minmaxyears, minmax)
}


###mean curve and spectrum of stack
mean_curves<-list()

for(j in 1:length(format_list)){
  df<-format_list[j]
  cluster_name<-names(df)
  df<-unlist(df,recursive=FALSE)
  pages_df<-bind_rows(df)
  pages_df$rwi<-as.numeric(pages_df$rwi)
  pages_df$year<-as.numeric(pages_df$year)
  pages_sum<-pages_df%>%
    group_by(year)%>%
    summarise(mean_rwi = mean(rwi, na.rm = T),
              n_sites = n())
  pages_sum<-as.data.frame(pages_sum)
  pages_sum$n_sites<-as.numeric(pages_sum$n_sites)
  mean_curve_spec<-as.data.frame(cbind(pages_sum$mean_rwi, pages_sum$year))
  colnames(mean_curve_spec)<-c("rwi","year")
  mean_curve_spec<-ts(mean_curve_spec$rwi, deltat = 1)
  mean_curve_spec<-SpecMTM(mean_curve_spec)
  mean_curve_spec$freq<-mean_curve_spec$freq[-(1:2)]
  mean_curve_spec$spec<-mean_curve_spec$spec[-(1:2)]
  mean_curve_spec$dof<-mean_curve_spec$dof[-(1:2)]
  smooth<-FilterSpec(mean_curve_spec, spans=c(3,5))
  smooth<-FilterSpecLog(smooth, df = 0.01)
  mean_curve<-as.data.frame(cbind(smooth$freq, smooth$spec, smooth$dof))
  colnames(mean_curve)<-c("freq","spectrum_of_stack","dofs_of_stack")
  mean_curve_final<-mean_curve %>%
    gather("type", "val", 2:3)
  colnames(mean_curve_final)<-c("freq",'type',"val")
  mean_curve_final<-list(mean_curve_final)
  names(mean_curve_final)<-cluster_name
  mean_curves<-c(mean_curves, mean_curve_final)
}

##################
theme1<-theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.background = element_rect(fill=alpha('white', 0)),
              legend.key.height = unit(0, 'lines'),
              legend.spacing.y = unit(0, 'cm'),
              text = element_text(size = 10),
              legend.position = c(0.7,0.2),
              legend.title = element_blank(),
              legend.direction="vertical",
              axis.text = element_text(color = "black"))



plot_list<-list()
curves_list<-list()
length_vec<-lapply(all_clusters, function(x) length(x))
plot_names<-names(mean_curves)                   

for(i in 1:length(mean_curves)) {
  cluster_name<-plot_names[[i]]
  length_clust<-paste("n =", as.character(length_vec)[[i]])
  mean_of_stack<-mean_curves[[i]]
  noiseandsignal<-spec_long_list[[i]]
  noiseandsignal<-noiseandsignal[,c("freq","val","type")]
  mean_of_stack<-mean_of_stack[,c("freq",'val',"type")]
  spec_df<-df_list[[i]]
  all_curves<-rbind(mean_of_stack, noiseandsignal)
  all_curves$freq<-round(all_curves$freq, digits = 8)
  curves_wide<-all_curves%>%
    group_by(as.factor(freq))%>%
    spread(type, val)%>%
    filter(!is.na(mean_of_all_spectra))
  curves_wide$`as.factor(freq)`<-NULL
  ############################################
  curves_wide$noise<-(curves_wide$mean_of_all_spectra - curves_wide$spectrum_of_stack)/(1-1/curves_wide$n.spec)
  curves_wide$signal<-curves_wide$mean_of_all_spectra - curves_wide$noise
  ##########################################
  signal<-list(curves_wide$freq, curves_wide$signal)
  names(signal)<-c("freq","spec")
  class(signal)<-'spectrum'
  signal<-FilterSpec(signal, spans = c(3,5))
  signal<-FilterSpecLog(signal, df = 0.1)
  curves_wide$signal<-signal$spec
  ##
  noise<-list(curves_wide$freq, curves_wide$noise)
  names(noise)<-c("freq","spec")
  class(noise)<-'spectrum'
  noise<-FilterSpec(noise, spans = c(3,5))
  noise<-FilterSpecLog(noise, df = 0.1)
  curves_wide$noise<-noise$spec
  ########################################calculate SNR
  curves_wide$ratio<-curves_wide$signal/curves_wide$noise
  #################################################
  p1<-ggplot()+
    geom_line(data = spec_df, aes(x=freq, y=spec, group = name, color = "detrending_type", linetype = "detrending_type"), size = 0.4)+
    geom_line(data = curves_wide, aes(x=freq, y=mean_of_all_spectra, color = "mean_of_all_spectra", linetype = "mean_of_all_spectra"), size = 1)+
    geom_line(data = curves_wide, aes(x=freq, y=signal, color = "signal", linetype = "signal"), size = 1)+
    geom_line(data = curves_wide, aes(x=freq, y=noise, color = "noise", linetype = "noise"), size = 1)+
    geom_line(data = curves_wide, aes(x = freq, y=spectrum_of_stack, color="spectrum_of_stack",linetype = "spectrum_of_stack"), size = 1)+
    scale_x_continuous(trans = c('log10'), limits = c(0.005,0.5),
                       breaks = c(0.01, 0.033, 0.1, 0.5),
                       labels = c("100","30","10","5"))+
    scale_y_continuous(trans = c('log10'))+
    scale_color_manual(name = " ", breaks = c("detrending_type",  "mean_of_all_spectra", "signal", "noise", "spectrum_of_stack"),
                       labels = c("Individual series", "Mean of all spectra","Signal","Noise", "Spectrum of stack"),
                       values = c("#666666","#666666","#e41a1c", "#a6761d","blue"))+
    scale_linetype_manual(name = " ", breaks = c("detrending_type",  "mean_of_all_spectra", "signal", "noise", "spectrum_of_stack"),
                          labels = c("Individual series", "Mean of all spectra","Signal","Noise","Spectrum of stack"),
                          values = c("dotted","solid","solid", "solid","solid"))+
    theme_bw()+
    annotate("text", x = 0.01, y = 0.1, label = length_clust)+
    theme1+
    theme(panel.border=element_blank(), axis.line=element_line())+
    coord_capped_cart(bottom="both", left="both", right = "both")+
    guides(fill = guide_legend(byrow = TRUE))+
    ylab(expression(PSD ^2~yr^-1))+
    xlab("Timescale (years)")+
    ggtitle(cluster_name)
  curves_wide$name<-cluster_name
  curves_wide<-as.data.frame(curves_wide)
  curves_wide<-list(curves_wide)
  names(curves_wide)<-cluster_name
  plot_list[[i]]<- p1
  curves_list<-c(curves_list, curves_wide)
}

names(plot_list)<-names(curves_list)

#graph each plot if desired
#plot_list[[1]]

#######approximate curves (uncorrected mean curve, signal, noise, and SNR) onto the same frequency axis for averageing across clusters
approx_list<-list()

for(i in 1:length(curves_list)){
  raw<-curves_list[[i]]
  name<-names(curves_list[i])
  raw<-list(raw$freq, raw$mean_of_all_spectra, raw$dofs_of_stack)
  names(raw)<-c("freq","spec","dof")
  class(raw)<-'spectrum'
  raw<-SpecApprox(raw)
  raw<-FilterSpec(raw, spans = c(3,5))
  raw<-FilterSpecLog(raw, df = 0.1)
  raw<-AddConfInterval(raw)
  
  spec<-curves_list[[i]]
  snr<-list(spec$freq, spec$ratio)
  names(snr)<-c("freq","spec")
  class(snr)<-'spectrum'
  snr<-SpecApprox(snr)
  
  signal<-list(spec$freq, spec$signal)
  names(signal)<-c("freq","spec")
  class(signal)<-'spectrum'
  #signal<-LogSmooth(signal, df.log = 0.1, removeFirst = 0, removeLast = 0)
  signal<-SpecApprox(signal)
  
  noise<-list(spec$freq, spec$noise)
  names(noise)<-c("freq","spec")
  class(noise)<-'spectrum'
  noise<-SpecApprox(noise)
  
  df<-as.data.frame(cbind(snr$freq, snr$spec, signal$spec, noise$spec, raw$spec, raw$lim.1, raw$lim.2))
  df$name<-name
  colnames(df)<-c("freq","snr","signal","noise","raw","raw_upperCI","raw_lowerCI", "name")
  df<-list(df)
  names(df)<-names(curves_list)[[i]]
  approx_list<-c(approx_list, df)
}

###bind together into one dataframe, calculate mean uncorrected, signal, noise, and SNR by proxy type
signal_curves <- do.call(rbind, approx_list)

proxy<-strsplit(signal_curves$name, "_")
proxy<-do.call(rbind, proxy)
proxy<-as.data.frame(proxy)
signal_curves$proxy<-proxy$V2
row.names(signal_curves)<-NULL

#in case any density curves need to be re-assigned to MXD
signal_curves$proxy[signal_curves$proxy=="delta Density"]<-"MXD"

#maintain only high frequencies for which all crns are represented
signal_curves<-signal_curves%>%
  filter(freq < 0.49)                   

#calculate per-cluster parameters (alpha and beta) for signal and noise estimates
signal_betas <- signal_curves %>%
  group_by(name)%>%
  do(mod = lm(log(signal)~log(freq), data = .))%>%
  mutate(signalBeta = summary(mod)$coefficients[2],
         signalAlpha = summary(mod)$coefficients[1]) %>%
  dplyr::select(-mod)

noise_betas <- signal_curves %>%
  group_by(name)%>%
  do(mod = lm(log(noise)~log(freq), data = .))%>%
  mutate(Slope = summary(mod)$coefficients[2],
         Alpha = summary(mod)$coefficients[1]) %>%
  dplyr::select(-mod)

signal_betas$noiseBeta<-noise_betas$Slope
signal_betas$noiseAlpha<-noise_betas$Alpha

#calculate cluster length
cluster_n<-lapply(all_clusters, length)
cluster_n<-as.numeric(cluster_n)


#by proxy average uncorrected, signal, noise, and SNR
signal_mean<-signal_curves %>%
  group_by(freq, proxy)%>%
  filter(!is.na(signal))%>%
  filter(!is.na(snr))%>%
  filter(!is.na(noise))%>%
  filter(!is.na(raw))%>%
  summarise(mean_raw_curve  = mean(raw),
            rawUpperCI = mean(raw_upperCI),
            rawLowerCI = mean(raw_lowerCI),
            mean_signal_curve = mean(signal),
            mean_ratio_curve = mean(snr), 
            mean_noise_curve = mean(noise),
            n = n())

mean_snr<-signal_curves %>%
  group_by(name)%>%
  filter(!is.na(signal))%>%
  filter(!is.na(snr))%>%
  filter(!is.na(noise))%>%
  filter(!is.na(raw))%>%
  summarise(mean_snr = mean(snr))


ggplot(signal_mean, aes(x=freq, y=mean_signal_curve, color = "signal"))+
  facet_grid(.~proxy)+
  geom_line()+
  geom_ribbon(data = signal_mean, aes(x = freq, ymin = rawLowerCI , ymax = rawUpperCI), color = NA, alpha = 0.5)+
  geom_line(data = signal_mean, aes(x = freq, y=mean_raw_curve, color = "uncorrected"))+
  geom_line(data = signal_mean, aes(x = freq, y=mean_noise_curve, color = "noise"))+
  geom_line(data = signal_mean, aes(x = freq, y=mean_ratio_curve, color = "SNR"))+
  scale_x_continuous(trans = c('log10'), limits = c(0.01,0.49), breaks = c(0.01, 0.033, 0.1, 0.5), labels = c("100","30","10","5"))+  
  scale_y_continuous(trans = c('log10'))


signal_mean<-as.data.frame(signal_mean)

############contain all cluster data in one dataframe, including all sites in cluster, the start, end year, length, 
############interseries correlation, mean_rbar, eps, and sample depth (only for n. america sites)

minmaxyears_df<-as.data.frame(minmaxyears)
minmaxyears_df<-t(minmaxyears_df)
rownames(minmaxyears_df)<-NULL

colnames(minmaxyears_df)<-c("start","end","interseries_cor","mean_ncor","mean_rbar","mean_eps") #min_year, max_year, cor, mean_ncor, mean_rbar, mean_eps
minmaxyears_df<-as.data.frame(minmaxyears_df)
minmaxyears_df$name<-unlist(plot_names)
clusters_sites<-lapply(all_clusters, names)
minmaxyears_df$clusters<-as.character(c(clusters_sites))
minmaxyears_df$length<-minmaxyears_df$end-minmaxyears_df$start

minmaxyears_df<-right_join(minmaxyears_df, mean_snr, by = "name")
minmaxyears_df<-right_join(minmaxyears_df, signal_betas, by = "name")

proxy<-strsplit(minmaxyears_df$name, "_")
proxy<-do.call(rbind, proxy)
proxy<-as.data.frame(proxy)
minmaxyears_df$proxy<-proxy$V2

n_clusters<-lapply(all_clusters, length)
minmaxyears_df$n_clusters<-as.numeric(n_clusters)

plot(minmaxyears_df$interseries_cor, minmaxyears_df$mean_snr)

#write.csv(minmaxyears_df, file = "all_params.csv")
#write.csv(signal_curves,"allcurves.csv")
#write.csv(signal_mean, "mean_curves.csv")
######################################################################
#remove all objects except mean curves and parameters
#set up bootstrapping for calculating confidence intervals on the signal, noise and SNR using synthetic timeseries

rm(list = setdiff(ls(), c("minmaxyears_df", "signal_curves", "signal_mean", "FilterSpec","FilterSpecLog", "SpecApprox")))

betas<-minmaxyears_df
all_curves<-signal_curves
mean_curves<-signal_mean

#################
#option to look at only the longest sites in the analysis
#longnames<-c("78_TRW","79_TRW","80_TRW","81_TRW","82_TRW","168_TRW","169_TRW","215_TRW")

#all_curves<-all_curves%>%
#  filter(name %in% longnames)

# betas<-betas%>%
#   filter(name %in% longnames)

# mean_curves<-all_curves %>%
#   group_by(freq, proxy)%>%
#   filter(!is.na(signal))%>%
#   filter(!is.na(snr))%>%
#   filter(!is.na(noise))%>%
#   filter(!is.na(raw))%>%
#   summarise(mean_raw_curve  = mean(raw),
#             rawUpperCI = mean(raw_upperCI),
#             rawLowerCI = mean(raw_lowerCI),
#             mean_signal_curve = mean(signal),
#             mean_ratio_curve = mean(snr), 
#             mean_noise_curve = mean(noise),
#             n = n())
###############################

#Note this takes a horribly long time, best run on 10 or 100 iterations during development
all_sims<-list()

for(x in 1:1000) {
  sims_all<-data.frame(matrix(ncol = 7, nrow  = 0))
  run<-as.character(x)
  colnames(sims_all)<-c("freq","noise","signal","ratio","name","proxy","run")
  for(i in 1:nrow(betas)){  
    flavor = betas$name[[i]]
    proxy = betas$proxy[[i]]
    length = betas$length[[i]]
    sig_beta<-betas$signalBeta[[i]]*-1
    noise_beta<-betas$noiseBeta[[i]]*-1
    sig_alpha<-exp(betas$signalAlpha[[i]])
    noise_alpha<-exp(betas$noiseAlpha[[i]])
    n_sites<-betas$n_clusters[[i]]
    #simulate a single spectra with the signal parameters
    simclim<-ts(SimPLS(N = length, beta = sig_beta, alpha = sig_alpha))
    noise_list<-list()
    noise_df<-data.frame(matrix(ncol = length, nrow = 0))
    for (j in seq_along(1:n_sites)){
      #for n sites in each cluster, simulate a series with real noise parameters and add to the signal spec 
      sim<-ts(SimPLS(N = length, beta = noise_beta, alpha = noise_alpha))
      simplus<-sim+simclim
      simspec<-SpecMTM(simplus)
      simspec<-FilterSpec(simspec, spans = c(3,5))
      simspec<-FilterSpecLog(simspec, df = 0.01)
      simspec<-SpecApprox(simspec)
      simspec<-list(simspec)
      #combine all site spectra into a stack
      noise_list<-c(noise_list, simspec)
      noise_df[nrow(noise_df) + 1,] = simplus
    }
    #calculate the SNRs with synthetic data
    colnames(noise_df)<-seq_along(1:length)
    noise_df<-ts(colMeans(noise_df))
    spectrum_of_stack<-SpecMTM(noise_df)
    spectrum_of_stack<-FilterSpec(spectrum_of_stack, spans = c(3,5))
    spectrum_of_stack<-FilterSpecLog(spectrum_of_stack, df = 0.01)
    spectrum_of_stack<-SpecApprox(spectrum_of_stack)
    freq<-spectrum_of_stack$freq
    noise_list<-lapply(noise_list, function(x) rbind(x$spec))
    noise_list<-do.call(rbind, noise_list)
    mean_of_spectra<-colMeans(noise_list)
    spec_ests<-cbind.data.frame(freq, spectrum_of_stack$spec, mean_of_spectra)
    colnames(spec_ests)<-c("freq","spectrum_of_stack","mean_of_spectra")
    #############################
    spec_ests$noise<-(spec_ests$mean_of_spectra - spec_ests$spectrum_of_stack)/(1-1/n_sites)
    spec_ests$signal<-spec_ests$mean_of_spectra - spec_ests$noise
    ##############################
    signal<-list(spec_ests$freq, spec_ests$signal)
    names(signal)<-c("freq","spec")
    class(signal)<-'spectrum'
    signal<-FilterSpec(signal, spans = c(3,5))
    signal<-FilterSpecLog(signal, df = 0.1)
    spec_ests$signal<-signal$spec
    ##additional smoothing
    noise<-list(spec_ests$freq, spec_ests$noise)
    names(noise)<-c("freq","spec")
    class(noise)<-'spectrum'
    noise<-FilterSpec(noise, spans = c(3,5))
    noise<-FilterSpecLog(noise, df = 0.1)
    spec_ests$noise<-noise$spec
    #################################
    spec_ests$ratio<-spec_ests$signal/spec_ests$noise
    spec_ests$name<-flavor
    spec_ests$proxy<-proxy
    spec_ests$spectrum_of_stack<-NULL
    spec_ests$mean_of_spectra<-NULL
    spec_ests$run<-run
    sims_all<-rbind.data.frame(sims_all, spec_ests)
  }
  sims_mean<-sims_all%>%
    group_by(freq, proxy)%>%
    summarise(signal = mean(signal),
              noise = mean(noise),
              ratio = mean(ratio))
  sims<-list(sims_mean)
  all_sims<-c(all_sims, sims)
}


sims_df<-do.call(rbind, all_sims)

sims_mean_df<-sims_df %>%
  group_by(freq, proxy)%>%
  summarise(signal_sim = mean(signal),
            noise_sim = mean(noise),
            ratio_sim = mean(ratio),
            signal_0.9 = quantile(signal, 0.9, na.rm = TRUE),
            signal_0.1 = quantile(signal, 0.1, na.rm = TRUE),
            noise_0.9 = quantile(noise, 0.9, na.rm = TRUE),
            noise_0.1 = quantile(noise, 0.1, na.rm = TRUE),
            snr_0.9 = quantile(ratio, 0.9, na.rm = TRUE),
            snr_0.1 = quantile(ratio, 0.1, na.rm = TRUE))

ggplot(sims_mean_df, aes(x=freq, y=signal_0.9))+
  geom_line()+
  facet_grid(.~proxy)+
  geom_line(data = sims_mean_df, aes(x=freq, y=signal_0.1))+
  geom_line(data = sims_mean_df, aes(x=freq, y = signal_sim), color = "red",linetype = "longdash")+
  scale_x_continuous(trans = c("log10"))+
  scale_y_continuous(trans = c("log10"))


sims_mean_df$snr_0.1<-sims_mean_df$snr_0.1/sims_mean_df$ratio_sim
sims_mean_df$snr_0.9<-sims_mean_df$snr_0.9/sims_mean_df$ratio_sim

sims_mean_df$signal_0.1<-sims_mean_df$signal_0.1/sims_mean_df$signal_sim
sims_mean_df$signal_0.9<-sims_mean_df$signal_0.9/sims_mean_df$signal_sim

sims_mean_df$noise_0.1<-sims_mean_df$noise_0.1/sims_mean_df$noise_sim
sims_mean_df$noise_0.9<-sims_mean_df$noise_0.9/sims_mean_df$noise_sim

sims_mean_df$freq<-as.character(as.numeric(sims_mean_df$freq))
mean_curves$freq<-as.character(as.numeric(mean_curves$freq))

bootstrapped_curves<-full_join(mean_curves, sims_mean_df, by = c("freq","proxy"))
bootstrapped_curves$freq<-as.numeric(as.character(bootstrapped_curves$freq))
##if opted for only long sites
#bootstrapped_curves_longsites<-full_join(mean_curves, sims_mean_df, by = c("freq","proxy"))


colnames(bootstrapped_curves)<-c("freq","proxy","mean_proxy_curve","proxyUpperCI","proxyLowerCI","signal_est","snr_est","noise_est","n_spec",
                        "signal_sim","noise_sim","snr_sim","signal_0.9","signal_0.1","noise_0.9",
                        "noise_0.1","snr_0.9","snr_0.1")


bootstrapped_curves<-bootstrapped_curves%>%
  filter(proxy == "TRW" | freq>=0.0133)%>%
  filter(proxy == "MXD" | freq>=0.011)

ggplot(bootstrapped_curves, aes(x=freq, y=signal_est, color = "signal"))+
  facet_grid(.~proxy)+
  geom_line()+
  geom_ribbon(data = bootstrapped_curves, aes(ymin =  signal_est*signal_0.1, ymax = signal_est*signal_0.9, fill = "signal"), color = NA, alpha = 0.5)+
  geom_line(data = bootstrapped_curves, aes(y=noise_est, color = "noise"))+
  geom_ribbon(data = bootstrapped_curves, aes(ymin =  noise_est*noise_0.1, ymax = noise_est*noise_0.9, fill = "noise"), color = NA, alpha = 0.5)+
  geom_line(data = bootstrapped_curves, aes(y=snr_est, color = "snr"))+
  geom_ribbon(data = bootstrapped_curves, aes(ymin =  snr_est*snr_0.1, ymax = snr_est*snr_0.9, fill = "snr"), color = NA, alpha = 0.5)+
  geom_line(data = bootstrapped_curves, aes(y=mean_proxy_curve, color = "proxy"))+
  geom_ribbon(data = bootstrapped_curves, aes(ymin =  proxyLowerCI, ymax = proxyUpperCI, fill = "proxy"), color = NA, alpha = 0.5)+
  scale_x_continuous(trans = c('log10'), limits = c(0.01,0.49), breaks = c(0.01, 0.033, 0.1, 0.5), labels = c("100","30","10","5"))+  
  scale_y_continuous(trans = c('log10'))

