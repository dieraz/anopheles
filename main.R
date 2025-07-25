library(stringr)       # String manipulation functions
library(raster)        # Raster data handling
library(sf)            # Simple Features for spatial vector data
library(sp)            # Spatial data classes and methods
library(dplyr)         # Data manipulation
library(RColorBrewer)  # Color palettes for plots
library(blockCV)       # Spatial cross-validation tools
library(dismo)         # Species distribution modeling
library(gbm)           # Generalized Boosted Regression Models
library(readxl)        # Excel file reading

# Load An. stephensi data and remove duplicates based on coordinates and time range
data_Ansteph2 = read.csv("data_Ansteph2.csv")
data_Ansteph2 = data_Ansteph2[!duplicated(data_Ansteph2[c("latitude","longitude","year_start","year_end")]),]

# Load background vector data (excluding An. stephensi and unknowns)
data_VM = read.csv("background_not_stephensi.csv")
background_V1 <- data_VM

# Crop global shapefile to study region
xmin = -18;xmax = 112; ymin = -100; ymax = 45
shape = crop(shapefile("world-administrative-boundaries.shp"), extent(xmin,xmax,ymin,ymax)) 

# Build convex hull around An. stephensi presence points to define study region
selectedNodes = unique(data_Ansteph2[,c("longitude","latitude")])
hull = chull(selectedNodes); hull = c(hull,hull[1])
p = Polygon(selectedNodes[hull,]); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
contourPolygons_df = st_as_sf(SpatialPolygonsDataFrame(sps, data.frame(ID=1:length(sps))))
center <- st_centroid(contourPolygons_df)
vertices <-  st_coordinates(contourPolygons_df)[,1:2]
vertices <-  st_as_sf(as.data.frame(vertices), coords = c("X", "Y"))
furthest <- max(st_distance(center, vertices))

# Buffer presence points to include potentially suitable background
buffer_coords = st_as_sf(data_Ansteph2, coords=c("longitude","latitude"))
buffer_pres1 = st_buffer(buffer_coords, furthest)
background_pot = background_V1
background_pot_sf = st_as_sf(background_pot, coords=c("longitude","latitude"))

# Retain background points within the broader buffer
detection = sapply(st_intersects(background_pot_sf, buffer_pres1),function(x){length(x)!=0})
indices0 = which(detection)
background_pot2 = background_pot[indices0,]
background_pot2_sf = st_as_sf(background_pot2, coords=c("longitude","latitude"))

# Remove background points too close to presences (to avoid contamination)
dist = 1  # distance in CRS units (likely degrees)
buffer_pres2 = st_buffer(buffer_coords, dist)
exclusion = sapply(st_intersects(background_pot2_sf, buffer_pres2),function(x){length(x)==0})
indices3 = sort(unique(which(exclusion)))

# Setup for BRT modeling
nruns = 10                      # Number of model runs
brt_model_scvs = list()        # Store models
AUCs = list()                  # Store AUC values
rasters = list()               # Store prediction rasters

# Define datasets and timeframes for environmental drivers
datasets = c("20crv3","20crv3-era5","20crv3-w5e5","gswp3-w5e5")
DATASETS = c("20CRv3","20CRv3-ERA5","20CRv3-W5E5","GSWP3-W5E5")
datyear = c("2015","2021","2019","2019")

# Load historical land use and vegetation data
croplands = brick("landuse/landuse-totals_histsoc_annual_1901_2021.nc",varname = "cropland_total")
pastures = brick("landuse/landuse-totals_histsoc_annual_1901_2021.nc",varname = "pastures")
rangelands = brick("landuse/LUH2_GCB2019/LUH2_GCB2019_range_remapcon_1901_2019.nc")
urbanAreas = brick("landuse/landuse-urbanareas_histsoc_annual_1901_2021.nc")
primaryForest = brick("landuse/LUH2_GCB2019/LUH2_GCB2019_primf_remapcon_1901_2019.nc")
secondaryForest = brick("landuse/LUH2_GCB2019/LUH2_GCB2019_secdf_remapcon_1901_2019.nc")


# Loop over climate-land use datasets to prepare predictors and run BRTs
for(dt in 1:length(datasets)){
  
  # Load time series of climate and vegetation covariates (monthly or annual)
  temperature = brick(paste0("ISIMIP3a/obsclim/",DATASETS[dt],"/",datasets[dt],"_obsclim_tas_1901_",datyear[dt],"_monmean.nc"))
  precipitation = brick(paste0("ISIMIP3a/obsclim/",DATASETS[dt],"/",datasets[dt],"_obsclim_pr_1901_",datyear[dt],"_monmean.nc"))
  relativehumidity = brick(paste0("ISIMIP3a/obsclim/",DATASETS[dt],"/",datasets[dt],"_obsclim_hurs_1901_",datyear[dt],"_monmean.nc"))
  soilmoist = brick(paste0("water_global/soilmoist/ISIMIP3a/JULES-ES-VN6P3/obsclim/jules-es-vn6p3_",datasets[dt],"_obsclim_histsoc_default_soilmoist_global_monthly_1901_",datyear[dt],"_yearmonmean.nc"), varname = "soilmoist")
  gpp = brick(paste0("biomes/gpp-total/ISIMIP3a/JULES-ES-VN6P3/obsclim/jules-es-vn6p3_",datasets[dt],"_obsclim_histsoc_default_gpp-total_global_annual_1901_",datyear[dt],".nc"))
  
  # Prepare storage objects for this dataset
  brt_model_scvs[[dt]] = list()
  AUCs[[dt]] = matrix(nrow=nruns, ncol = 1)
  colnames(AUCs[[dt]]) = "SCV_AUC"
  
  numreps = 3  # Number of background replicates (3x presences)
  
  # Loop through multiple repetitions (different background draws)
  for(j in 1:nruns){
    selected = sample(length(indices3), numreps * nrow(data_Ansteph2), replace=FALSE)
    indices4 = indices3[selected]
    
    background_final = background_V1[rownames(background_V1) %in% indices4, ]
    row.names(background_final) = 1:nrow(background_final)
    
    # Combine presence and background points, assign response variable
    data = rbind(
      cbind(data_Ansteph2[,c("country","longitude","latitude","year_start","year_end")], rep(1, nrow(data_Ansteph2))),
      cbind(background_final[,c("country","longitude","latitude","year_start","year_end")], rep(0, nrow(background_final)))
    )
    colnames(data)[6] = "response"
    
    # Create empty data frame for covariates
    data_for_brt = as.data.frame(matrix(NA, nrow=nrow(data), ncol=17))
    data_for_brt[,1:6] = data[,1:6]  # Copy metadata
    colnames(data_for_brt) = c("country","longitude","latitude","year_start","year_end","response",
                               "temperature","precipitation","relative_humidity",
                               "croplands","pastures","rangelands","urbanAreas",
                               "primaryForest","secondaryForest","soilmoist","gpp")
    
    # Optional filter to exclude post-2016 observations for 20CRv3
    if(dt == 1){
      data_for_brt = data_for_brt[data_for_brt$year_end < 2016 & data_for_brt$year_start < 2016, ]
    }
    
    # Loop through all records and extract mean climate and land-use variables
    for(k in 1:nrow(data_for_brt)){
      year_start = data_for_brt[k, "year_start"]
      year_end = data_for_brt[k, "year_end"]
      month_start = (year_start - 1901)*12 + 1
      month_end = (year_end - 1901)*12 + 12
      
      temperature_temp = mean(temperature[[month_start:month_end]]) - 273.15  # Kelvin to Celsius
      precipitation_temp = mean(precipitation[[month_start:month_end]]) * 60*60*24  # mm/s to mm/day
      relhumidity_temp = mean(relativehumidity[[month_start:month_end]])
      
      year_start_temp = year_start - 1900
      year_end_temp = year_end - 1900
      
      croplands_temp = mean(croplands[[year_start_temp:year_end_temp]])
      pastures_temp = mean(pastures[[year_start_temp:year_end_temp]])
      rangelands_temp = mean(rangelands[[year_start_temp:year_end_temp]])
      urbanAreas_temp = mean(urbanAreas[[year_start_temp:year_end_temp]])
      primaryForest_temp = mean(primaryForest[[year_start_temp:year_end_temp]])
      secondaryForest_temp = mean(secondaryForest[[year_start_temp:year_end_temp]])
      soilmoist_temp = mean(soilmoist[[year_start_temp:year_end_temp]])
      gpp_temp = mean(gpp[[year_start_temp:year_end_temp]]) * 60*60*24  # convert to daily equivalent
      
      # Extract values from rasters at each point
      lonlat = data_for_brt[k, c("longitude", "latitude")]
      data_for_brt[k, "temperature"] = raster::extract(temperature_temp, lonlat)
      data_for_brt[k, "precipitation"] = raster::extract(precipitation_temp, lonlat)
      data_for_brt[k, "relative_humidity"] = raster::extract(relhumidity_temp, lonlat)
      data_for_brt[k, "croplands"] = raster::extract(croplands_temp, lonlat)
      data_for_brt[k, "pastures"] = raster::extract(pastures_temp, lonlat)
      data_for_brt[k, "rangelands"] = raster::extract(rangelands_temp, lonlat)
      data_for_brt[k, "urbanAreas"] = raster::extract(urbanAreas_temp, lonlat)
      data_for_brt[k, "primaryForest"] = raster::extract(primaryForest_temp, lonlat)
      data_for_brt[k, "secondaryForest"] = raster::extract(secondaryForest_temp, lonlat)
      data_for_brt[k, "soilmoist"] = raster::extract(soilmoist_temp, lonlat)
      data_for_brt[k, "gpp"] = raster::extract(gpp_temp, lonlat)
    }
    
    # Remove rows with missing values before model fitting
    data = na.omit(data_for_brt)
    
    # Define BRT settings
    gbm.x = 7:17  # columns with environmental predictors
    gbm.y = "response"
    n.folds = 5
    
    # Spatial cross-validation: partition study area into spatial blocks
    spdf = SpatialPointsDataFrame(data[,c("longitude","latitude")], data, proj4string = crs(shape))
    myblocks = cv_spatial(spdf, column="response", k=n.folds, raster=shape, size=c(1200,1200)*1000, selection="random")
    fold.vector = myblocks$folds_ids
    
    # Train the BRT model with spatial cross-validation
    brt_model_scvs[[dt]][[j]] = gbm.step(
      data = data,
      gbm.x = gbm.x,
      gbm.y = gbm.y,
      fold.vector = fold.vector,
      tree.complexity = 5,
      learning.rate = 0.005,
      bag.fraction = 0.80,
      family = "bernoulli",
      n.trees = 100,
      step.size = 10,
      max.trees = 10000,
      tolerance.method = "auto",
      tolerance = 0.001,
      verbose = TRUE
    )
    
    # Store cross-validation AUC
    AUCs[[dt]][j, "SCV_AUC"] = brt_model_scvs[[dt]][[j]]$cv.statistics$discrimination.mean
  }
}


# -----------  Predictions in the past -----------

# This block initializes parameters for generating predictions across multiple time periods
# and comparing counterfactual (no climate change) vs observed scenarios
models = c("obsclim","counterclim")  # Observed and counterfactual climate models
datasets = c("20crv3","20crv3-era5","20crv3-w5e5","gswp3-w5e5")
DATASETS = c("20CRv3","20CRv3-ERA5","20CRv3-W5E5","GSWP3-W5E5")
datyear = c("2015","2021","2019","2019")  # Latest year in each dataset

# Define historical periods to generate predictions for
years_start = c(1901,1920,1940,1960,1980,2000)
years_end = c(1919,1939,1959,1979,1999,2019)
years_end2 = c(2015,2019,2019,2019)  # actual dataset-dependent end years

# Load land use layers for the prediction phase (if not already loaded)
croplands = brick("landuse/landuse-totals_histsoc_annual_1901_2021.nc",varname = "cropland_total")
pastures = brick("landuse/landuse-totals_histsoc_annual_1901_2021.nc",varname = "pastures")
rangelands = brick("landuse/LUH2_GCB2019/LUH2_GCB2019_range_remapcon_1901_2019.nc")
urbanAreas = brick("landuse/landuse-urbanareas_histsoc_annual_1901_2021.nc")
primaryForest = brick("landuse/LUH2_GCB2019/LUH2_GCB2019_primf_remapcon_1901_2019.nc")
secondaryForest = brick("landuse/LUH2_GCB2019/LUH2_GCB2019_secdf_remapcon_1901_2019.nc")

# Initialize result containers
rasters_obsclim = list()
rasters_counterclim = list()

# Loop over datasets and models (obsclim vs counterclim)
for(dt in 3:length(datasets)){
  for(m in 1:length(models)){
    
    # Load climatic drivers for selected model
    temperature = brick(paste0("ISIMIP3a/",models[m],"/",DATASETS[dt],"/",datasets[dt],"_",models[m],"_tas_1901_",datyear[dt],"_monmean.nc"))
    precipitation = brick(paste0("ISIMIP3a/",models[m],"/",DATASETS[dt],"/",datasets[dt],"_",models[m],"_pr_1901_",datyear[dt],"_monmean.nc"))
    relativehumidity = brick(paste0("ISIMIP3a/",models[m],"/",DATASETS[dt],"/",datasets[dt],"_",models[m],"_hurs_1901_",datyear[dt],"_monmean.nc"))
    soilmoist = brick(paste0("water_global/soilmoist/ISIMIP3a/JULES-ES-VN6P3/",models[m],"/jules-es-vn6p3_",datasets[dt],"_",models[m],"_histsoc_default_soilmoist_global_monthly_1901_",datyear[dt],"_yearmonmean.nc"), varname = "soilmoist")
    gpp = brick(paste0("biomes/gpp-total/ISIMIP3a/JULES-ES-VN6P3/",models[m],"/jules-es-vn6p3_",datasets[dt],"_",models[m],"_histsoc_default_gpp-total_global_annual_1901_",datyear[dt],".nc"))
    
    # Initialize per-model list to hold rasters
    if(m == 1){rasters_obsclim[[dt]] = list()}
    if(m == 2){rasters_counterclim[[dt]] = list()}
    
    # Loop over defined historical periods
    for(y in 1:length(years_start)){
      
      if(m == 1){rasters_obsclim[[dt]][[y]] = list()}
      if(m == 2){rasters_counterclim[[dt]][[y]] = list()}
      
      # Loop over BRT replicates
      for(j in 1:length(brt_model_scvs[[dt]])){
        
        model = brt_model_scvs[[dt]][[j]]
        n.trees = model$gbm.call$best.trees
        
        # Define the start/end dates for this time window
        year_start = years_start[y]
        year_end = ifelse(y == 6, years_end2[dt], years_end[y])
        
        # Temporal slicing of predictors (climate and land use)
        month_start = (year_start - 1901)*12+1
        month_end = (year_end - 1901)*12+12
        temperature_temp = mean(temperature[[month_start:month_end]]) - 273.15
        precipitation_temp = mean(precipitation[[month_start:month_end]]) * 60*60*24
        relhumidity_temp = mean(relativehumidity[[month_start:month_end]])
        
        year_start_temp = year_start - 1900
        year_end_temp = year_end - 1900
        croplands_temp = mean(croplands[[year_start_temp:year_end_temp]])
        pastures_temp = mean(pastures[[year_start_temp:year_end_temp]])
        rangelands_temp = mean(rangelands[[year_start_temp:year_end_temp]])
        urbanAreas_temp = mean(urbanAreas[[year_start_temp:year_end_temp]])
        primaryForest_temp = mean(primaryForest[[year_start_temp:year_end_temp]])
        secondaryForest_temp = mean(secondaryForest[[year_start_temp:year_end_temp]])
        soilmoist_temp = mean(soilmoist[[year_start_temp:year_end_temp]])
        gpp_temp = mean(gpp[[year_start_temp:year_end_temp]]) * 60*60*24
        
        # Stack covariates into one raster object
        envVariables = list(temperature_temp, precipitation_temp, relhumidity_temp,
                            croplands_temp, pastures_temp, rangelands_temp, urbanAreas_temp,
                            primaryForest_temp, secondaryForest_temp, soilmoist_temp, gpp_temp)
        
        # Crop/mask all rasters to study region
        for (i in 1:length(envVariables)) envVariables[[i]] = crop(envVariables[[i]], shape, snap="out")
        for (i in 1:length(envVariables)) envVariables[[i]] = mask(envVariables[[i]], shape)
        
        # Align NA values across layers
        NAcells = unique(unlist(lapply(envVariables, function(x) which(is.na(x[])))))
        for (i in 1:length(envVariables)) envVariables[[i]][NAcells] = NA
        
        # Stack covariates and make prediction
        rasters_stack = stack(envVariables)
        df = as.data.frame(rasters_stack)
        not_NA = which(!is.na(rowMeans(df)))
        newdata = df[not_NA,]
        colnames(newdata) = colnames(data_for_brt)[7:17]
        
        prediction = predict.gbm(model, newdata, n.trees, type = "response")
        
        # Write predicted values to raster
        rast = envVariables[[1]]
        rast[which(!is.na(rast[]))] = prediction
        
        if(m == 1){rasters_obsclim[[dt]][[y]][[j]] = rast}
        if(m == 2){rasters_counterclim[[dt]][[y]][[j]] = rast}
      }
    }
  }
}


# -----------  Plots  -----------

# Generate PDF maps of predictions for past periods (obsclim and counterclim scenarios)
legend1 = raster(as.matrix(c(0,1)))  # Legend raster for predicted suitability values

# Loop over datasets to create comparative plots
for(dt in 1:length(datasets)){
  
  pdf(paste0("figures/past_",datasets[dt],".pdf"), width=12, height=4)
  par(mfrow=c(2,6), oma=c(0,0,1.5,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
  
  # Plot observed climate predictions (top row)
  for(i in 1:length(rasters_obsclim[[dt]])){
    raster2plot = mean(stack(rasters_obsclim[[dt]][[i]]))
    plot(raster2plot, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131][1:(max(raster2plot[],na.rm=T)*100)],
         ann=F, legend=F, axes=F, box=F)
    mtext(paste0(years_start[i],"-",years_end[i]), side=3, line=-0.7, cex=0.65, col="gray30")
    if(i == 6){
      plot(legend1, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131], legend.only=T, add=T,
           legend.width=0.5, legend.shrink=0.3,
           smallplot=c(0.10,0.80,0.20,0.23), adj=3,
           axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30",
                          tck=-0.6, col.axis="gray30", line=0, mgp=c(0,0.0,0),
                          at=seq(0,1,0.25), labels=c("0","0.25","0.5","0.75","1")))
    }
  }
  
  # Plot counterfactual climate predictions (bottom row)
  for(i in 1:length(rasters_counterclim[[dt]])){
    raster2plot = mean(stack(rasters_counterclim[[dt]][[i]]))
    plot(raster2plot, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131][1:(max(raster2plot[],na.rm=T)*100)],
         ann=F, legend=F, axes=F, box=F)
    mtext(paste0(years_start[i],"-",years_end[i]), side=3, line=-0.7, cex=0.65, col="gray30")
  }
  
  dev.off()
}


# -----------  Climate Attribution Visualization -----------

# Create a figure comparing early vs late period and their climate attribution difference
legend1 = raster(as.matrix(c(0,1)))  # Legend for predicted suitability
legend2 = raster(as.matrix(c(-0.5,0.5)))  # Legend for difference maps

pdf("figures/figure1.pdf", width=6, height=3)
par(mfrow=c(2,3), oma=c(0,0,1.5,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")

# 1. Early obsclim map (e.g., 1901–1919)
raster2plot = mean(stack(rasters_obsclim[[4]][[1]]))
plot(raster2plot, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131][1:(max(raster2plot[],na.rm=T)*100)],
     ann=F, legend=F, axes=F, box=F)
mtext(paste0(years_start[1],"-",years_end[1]), side=3, line=-0.3, cex=0.65, col="gray30")

# 2. Late obsclim map (e.g., 2000–2019)
raster2plot = mean(stack(rasters_obsclim[[4]][[6]]))
plot(raster2plot, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131][1:(max(raster2plot[],na.rm=T)*100)],
     ann=F, legend=F, axes=F, box=F)
mtext(paste0(years_start[6],"-",years_end[6]), side=3, line=-0.3, cex=0.65, col="gray30")

# 3. Climate effect map: difference between obsclim and counterclim
raster2plot = mean(stack(rasters_obsclim[[4]][[6]])) - mean(stack(rasters_counterclim[[4]][[6]]))
# Normalize the difference map to 0-1 scale
min_val <- min(raster2plot[], na.rm = TRUE)
max_val <- max(raster2plot[], na.rm = TRUE)
raster2plot_normalized <- (raster2plot - min_val) / (max_val - min_val)
plot(raster2plot_normalized, col = rev(colorRampPalette(brewer.pal(11, "BrBG"))(131)), ann = FALSE, legend = FALSE, axes = FALSE, box = FALSE)
mtext("difference", side=3, line=-0.3, cex=0.65, col="gray30")

# 4. Early counterclim map
raster2plot = mean(stack(rasters_counterclim[[4]][[1]]))
plot(raster2plot, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131][1:(max(raster2plot[],na.rm=T)*100)],
     ann=F, legend=F, axes=F, box=F)
mtext(paste0(years_start[1],"-",years_end[1]), side=3, line=-0.3, cex=0.65, col="gray30")

# 5. Late counterclim map
raster2plot = mean(stack(rasters_counterclim[[4]][[6]]))
plot(raster2plot, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131][1:(max(raster2plot[],na.rm=T)*100)],
     ann=F, legend=F, axes=F, box=F)
mtext(paste0(years_start[6],"-",years_end[6]), side=3, line=-0.3, cex=0.65, col="gray30")

# 6. Add legends for prediction and difference
plot(legend1, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131], legend.only=T, add=F, legend.width=0.5,
     legend.shrink=0.3, smallplot=c(0.10,0.13,0.20,0.9), adj=1,
     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col.axis="gray30",
                    line=0, mgp=c(0,0.0,0), at=seq(0,1,0.25), labels=c("0","0.25","0.5","0.75","1")))

plot(legend2, col=rev(colorRampPalette(brewer.pal(11,"BrBG"))(131))[21:131], legend.only=T, add=F, legend.width=0.5,
     legend.shrink=0.3, smallplot=c(0.20,0.23,0.20,0.9), adj=1,
     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col.axis="gray30",
                    line=0, mgp=c(0,0.0,0), at=seq(-0.5,0.5,0.25), labels=c("-0.5","-0.25","0","0.25","0.5")))

dev.off()


# -----------  Relative Influence Analysis -----------

# Define variable names used in the models
envVariableNames = c("temperature","precipitation","relative_humidity",
                     "croplands","pastures","rangelands","urbanAreas",
                     "primaryForest","secondaryForest","soilmoist","gpp")

# Loop over datasets to compute relative influence of each predictor
for(dt in 1:length(datasets)){
  relativeInfluences = matrix(0, nrow=length(envVariableNames), ncol=3)  # Mean and CI bounds
  relativeInfluences_all = matrix(0, nrow=length(envVariableNames), ncol=nruns)
  
  brt_model_scv = brt_model_scvs[[dt]]  # Get replicate models for one dataset
  
  # Accumulate relative influence per model and variable
  for (j in 1:length(brt_model_scv)){
    for (k in 1:length(envVariableNames)){
      relinf = summary(brt_model_scv[[j]])[envVariableNames[k], "rel.inf"]
      relativeInfluences[k,1] = relativeInfluences[k,1] + relinf
      relativeInfluences_all[k,j] = relinf
    }
  }
  
  # Normalize mean influence and compute confidence intervals
  row.names(relativeInfluences) = envVariableNames
  relativeInfluences[,1] = relativeInfluences[,1]/length(brt_model_scv)
  relativeInfluences[,2] = apply(relativeInfluences_all, 1, function(x) quantile(x, 0.025, na.rm = TRUE))
  relativeInfluences[,3] = apply(relativeInfluences_all, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
  
  # Extract median, min, and max values of variables at presence locations
  data_curv = data[data$response == 1, ]
  envVariableValues = matrix(nrow=3, ncol=length(envVariableNames))
  row.names(envVariableValues) = c("median","minV","maxV")
  colnames(envVariableValues) = envVariableNames
  for (j in 1:length(envVariableNames)){
    envVariableValues[,j] = c(median(data_curv[,envVariableNames[j]], na.rm=T),
                              min(data_curv[,envVariableNames[j]], na.rm=T),
                              max(data_curv[,envVariableNames[j]], na.rm=T))
  }
  
  # Plot partial dependence curves
  pdf(paste0("figures/relative_influences_",datasets[dt],".pdf"), width=8, height=6)
  par(mfrow=c(3,4), oma=c(1.3,1.5,1,0.5), mar=c(2.5,1,0.5,1), lwd=0.2, col="gray30")
  
  for (i in 1:length(envVariableNames)){
    df = data.frame(matrix(nrow=100, ncol=length(envVariableNames)))
    colnames(df) = envVariableNames
    for (k in 1:length(envVariableNames)){
      if (i == k){
        df[,envVariableNames[k]] = seq(envVariableValues["minV",k], envVariableValues["maxV",k], length.out=100)
      } else {
        df[,envVariableNames[k]] = rep(envVariableValues["median",k], 100)
      }
    }
    
    # Predict using each BRT model
    predictions_list = list()
    for (j in 1:length(brt_model_scv)){
      n.trees = brt_model_scv[[j]]$gbm.call$best.trees
      pred = predict.gbm(brt_model_scv[[j]], newdata=df, n.trees=n.trees, type="response")
      predictions_list[[j]] = pred
    }
    
    # Determine plotting range
    minX = min(df[,envVariableNames[i]])
    maxX = max(df[,envVariableNames[i]])
    minY = min(unlist(predictions_list))
    maxY = max(unlist(predictions_list))
    
    # Plot curves
    for (l in 1:length(predictions_list)){
      if (l == 1){
        plot(df[,envVariableNames[i]], predictions_list[[l]], type="l", col="red", lwd=0.2,
             ann=F, axes=F, xlim=c(minX, maxX), ylim=c(minY, maxY))
      } else {
        lines(df[,envVariableNames[i]], predictions_list[[l]], col="red", lwd=0.2)
      }
    }
    
    # Format axes and titles
    envVariableNames1 = c("temperature","precipitation","relative humidity",
                          "croplands","pastures","rangelands","urban areas",
                          "primary forest","secondary forest","soil moisture",
                          "gross primary productivity")
    axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.07,0))
    axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.2,0))
    title(ylab="predicted values", cex.lab=0.9, mgp=c(1.3,0,0), col.lab="gray30")
    title(xlab=paste0(envVariableNames1[i], " (", round(relativeInfluences[i,1],1), "% [",
                      round(relativeInfluences[i,2],1), "-", round(relativeInfluences[i,3],1), "])"),
          cex.lab=0.9, mgp=c(1.3,0,0), col.lab="gray30")
    box(lwd=0.2, col="gray30")
  }
  dev.off()
}


# -----------  Sørensen Index Evaluation -----------

# Calculate the prevalence-pseudoabsence-calibrated Sørensen index (SIppc) to evaluate model performance
# Based on: Leroi et al. (2018, J. Biogeography) and Li & Guo (2013, Ecography)

# Initialize storage for results
tabs_list1 = list()  # stores Sørensen curves for each dataset
SIppcs = matrix(nrow=nruns, ncol=length(datasets))
colnames(SIppcs) = datasets
thresholds = matrix(nrow=nruns, ncol=length(datasets))
colnames(thresholds) = datasets

# Loop over datasets
for (i in 1:length(datasets)) {
  tabs_list2 = list()
  brt_model_scv = brt_model_scvs[[i]]
  
  # Evaluate each replicate model
  for (j in 1:length(brt_model_scv)) {
    model = brt_model_scv[[j]]
    df = model$gbm.call$dataframe
    responses = df$response
    predictors = df[,2:ncol(df)]
    n.trees = model$gbm.call$best.trees
    
    # Generate predictions for internal data
    prediction = predict.gbm(model, predictors, n.trees, type = "response")
    P = sum(responses == 1)
    A = sum(responses == 0)
    prev = P / (P + A)
    x = (P / A) * ((1 - prev) / prev)  # correction factor for prevalence
    
    # Search for threshold that maximizes SIppc
    best_sorensen = 0
    best_threshold = 0
    tmp = matrix(nrow=101, ncol=2)
    tmp[,1] = seq(0, 1, 0.01)
    
    for (threshold in seq(0, 1, 0.01)) {
      TP = sum((responses == 1) & (prediction >= threshold))
      FN = sum((responses == 1) & (prediction < threshold))
      FP_pa = sum((responses == 0) & (prediction >= threshold))
      
      sorensen = (2 * TP) / ((2 * TP) + (x * FP_pa) + FN)
      tmp[which(tmp[,1] == threshold), 2] = sorensen
      
      if (sorensen > best_sorensen) {
        best_sorensen = sorensen
        best_threshold = threshold
      }
    }
    
    tabs_list2[[j]] = tmp
    SIppcs[j, i] = best_sorensen
    thresholds[j, i] = best_threshold
  }
  
  tabs_list1[[i]] = tabs_list2
}

# Plot the Sørensen index profiles
pdf("malaria/figures/SI_ppcs.pdf", width=8, height=1.8)
par(mfrow=c(1,4), oma=c(0,0.3,0,0), mar=c(2.5,2.5,0.5,0.5), lwd=0.4, col="gray30")

for (i in 1:length(datasets)) {
  plot(tabs_list1[[i]][[1]], col=NA, ann=F, axes=F, xlim=c(0,1), ylim=c(0,1))
  for (j in 1:length(tabs_list1[[i]])) {
    lines(tabs_list1[[i]][[j]], lwd=0.3, col="gray80")
  }
  axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.025, col.axis="gray30", mgp=c(0,0.14,0))
  axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.025, col.axis="gray30", mgp=c(0,0.35,0))
  if (i == 1) {
    title(ylab=expression("SI"["ppc"]), cex.lab=0.9, mgp=c(1.3,0,0), col.lab="gray30")
  }
  title(xlab="threshold", cex.lab=0.9, mgp=c(1.1,0,0), col.lab="gray30")
  box(lwd=0.2, col="gray30")
  mtext(datasets[i], side=3, line=-1.3, at=0.995, cex=0.55, col="gray30", adj=1)
}

dev.off()

