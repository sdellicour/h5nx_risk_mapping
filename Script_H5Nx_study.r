library(ade4)
library(ape)
library(blockCV)
library(diagram)
library(dismo)
library(exactextractr)
library(fields)
library(gbm)
library(geosphere)
library(ggplot2)
library(HDInterval)
library(lubridate)
library(ncdf4)
library(ncf)
library(picante)
library(pgirmess)
library(phytools)
library(pROC)
library(RColorBrewer)
library(raster)
library(rgdal)
library(rgeos)
library(seegSDM)
library(seqinr)
library(sf)
library(sp)
library(vioplot)

# 1. Preparation of the environmental rasters
# 2. Loading and expecting the occurrence records
# 3. Preparing the BRT dataframe for each dataset
# 4. All boosted regression trees (BRT) analyses
# 5. Comparing all the AUC support and SSB values
# 6. Computing the relative influence (RI) values
# 7. Comparing the predictive capacities (with AUC)
# 8. Comparing the response curves and RI values
# 9. Projections of H5 ecological suitabilities

savingPlots = FALSE; showingPlots = FALSE; writingFiles = FALSE

# 1. Preparation of the environmental rasters

	# 1.1. Generating the mask rasters

if (!file.exists("Environmental_data/Prepared_rasters/Mask_raster_08.tif"))
	{
		if (!file.exists("Environmental_data/Prepared_rasters/Mask_raster_008.tif"))
			{
				buffer1 = raster("Environmental_data/Original_rasters/Day_LST_annual_mean.tif")
				buffer1[!is.na(buffer1[])] = 0
				writeRaster(buffer1, "Environmental_data/Prepared_rasters/Mask_raster_008.tif")
			}
		buffer1 = raster("Environmental_data/Prepared_rasters/Mask_raster_008.tif")
		buffer2 = buffer1; buffer2[!is.na(buffer2[])] = 1; buffer2[is.na(buffer2[])] = 0
		buffer3 = raster::aggregate(buffer2, 10, fun=mean)
		buffer4 = buffer3; buffer4[buffer4[]<0.5] = NA; buffer4[buffer4[]>=0.5] = 0
		lakes = shapefile("Environmental_data/NaturalEarth_files/Natural_Earth_lakes.shp")
		for (i in 1:length(lakes@polygons)) # to discard the lakes from the mask
			{
				for (j in 1:length(lakes@polygons[[i]]@Polygons))
					{
						pol = lakes@polygons[[i]]@Polygons[[j]]
						p = Polygon(pol@coords); ps = Polygons(list(p),1)
						sps = SpatialPolygons(list(ps)); pol = sf::st_as_sfc(sps)
						extraction = exact_extract(buffer4, pol, include_cell=T)[[1]]
						extraction = extraction[which(extraction[,"coverage_fraction"]>0.5),]
						buffer4[extraction[,"cell"]] = NA	
					}
			}
		caspianSea = shapefile("Environmental_data/NaturalEarth_files/Caspian_sea_pol.shp")	
		for (i in 1:length(caspianSea@polygons)) # to discard the Caspian Sea from the mask
			{
				for (j in 1:length(caspianSea@polygons[[i]]@Polygons))
					{
						pol = caspianSea@polygons[[i]]@Polygons[[j]]
						p = Polygon(pol@coords); ps = Polygons(list(p),1)
						sps = SpatialPolygons(list(ps)); pol = sf::st_as_sfc(sps)
						extraction = exact_extract(buffer4, pol, include_cell=T)[[1]]
						extraction = extraction[which(extraction[,"coverage_fraction"]>0.5),]
						buffer4[extraction[,"cell"]] = NA	
					}
			}
		writeRaster(buffer4, "Environmental_data/Prepared_rasters/Mask_raster_08.tif") 
	}

	# 1.2. Uniformising the extent and resolution of all rasters

mask = raster("Environmental_data/Prepared_rasters/Mask_raster_08.tif")
data_sources = read.csv("Environmental_data/Data_sources.csv", head=T, sep=";")
envVariableNames = unique(data_sources[,"variable_name_2"]); envVariableNames = envVariableNames[envVariableNames!=""]
envVariableNames = envVariableNames[!envVariableNames%in%c("vaccination_in_china")]
for (i in 1:length(envVariableNames))
	{
		if (!file.exists(paste0("Environmental_data/Prepared_rasters/",envVariableNames[i],".tif")))
			{
				envVariableName = envVariableNames[i]; substr(envVariableName,1,1) = toupper(substr(envVariableName,1,1))				
				rast = crop(raster(paste0("Environmental_data/Original_rasters/",envVariableName,".tif")), extent(-180,180,-56,90))
				if (round(res(rast)[1],5) == 0.00833)
					{
						fun = data_sources[which(data_sources[,"variable_name_2"]==envVariableNames[i])[1],"aggregation_function"]
						if (fun == "sum")
							{
								rast[is.na(rast[])] = 0
								rast = raster::aggregate(rast, 10, fun=sum)
							}
						if (fun == "mean1")
							{
								rast[is.na(rast[])] = 0
								rast = raster::aggregate(rast, 10, fun=mean)
							}
						if (fun == "mean2")
							{
								rast = raster::aggregate(rast, 10, fun=mean, na.rm=T)
							}
					}
				rast[is.na(mask)] = NA
				writeRaster(rast, paste0("Environmental_data/Prepared_rasters/",envVariableName,".tif"))
			}
	}

	# 1.3. Plotting all the environmental variables

allEnvVariableNames = unique(data_sources[,"variable_name_2"]); allEnvVariableNames = allEnvVariableNames[allEnvVariableNames!=""]
allEnvVariableNames = allEnvVariableNames[!allEnvVariableNames%in%c("vaccination_in_china")]
envVariableNamesPlot = allEnvVariableNames[!allEnvVariableNames%in%c("open_water_areas","urban_and_built_up_areas")]
allEnvVariables = list(); envVariablesPlot = list()
for (i in 1:length(allEnvVariableNames))
	{
		envVariableName = allEnvVariableNames[i]; substr(envVariableName,1,1) = toupper(substr(envVariableName,1,1))				
		allEnvVariables[[i]] = raster(paste0("Environmental_data/Prepared_rasters/",envVariableName,".tif"))
		if (grepl("population_density",envVariableName))
			{
				allEnvVariables[[i]][] = log(allEnvVariables[[i]][]+1)
				allEnvVariableNames[i] = paste0(allEnvVariableNames[i],"_(log)")
			}
	}
for (i in 1:length(envVariableNamesPlot))
	{
		envVariableName = envVariableNamesPlot[i]; substr(envVariableName,1,1) = toupper(substr(envVariableName,1,1))				
		envVariablesPlot[[i]] = raster(paste0("Environmental_data/Prepared_rasters/",envVariableName,".tif"))
		if (grepl("population_density",envVariableName))
			{
				envVariablesPlot[[i]][] = log(envVariablesPlot[[i]][]+1)
				envVariableNamesPlot[i] = paste0(envVariableNamesPlot[i],"_(log)")
			}
	}
if (savingPlots)
	{
		coast_lines = gSimplify(crop(shapefile("Environmental_data/NaturalEarth_files/Coastline_borders.shp"), extent(-180,180,-56,90)), 0.01)
		colourScale = rev(colorRampPalette(brewer.pal(9,"RdYlBu"))(120)[10:110])
		i1s = c(1,37); i2s = c(36,length(envVariablesPlot))
		for (h in 1:length(i1s))
			{		
				pdf(paste0("Env_variables_",h,".pdf"), width=14, height=18) # dev.new(width=14, height=18)
				par(mfrow=c(9,4), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30"); options(scipen=9)
				for (i in i1s[h]:i2s[h])
					{
						plot(envVariablesPlot[[i]], col=colourScale, ann=F, axes=F, box=F, legend=F)
						plot(coast_lines, lwd=0.1, col="gray30", add=T)
						envVariableName = gsub("_2010_\\(","_\\(2010,_",envVariableNamesPlot[i])
						envVariableName = gsub("_2015_\\(","_\\(2015,_",envVariableName)
						envVariableName = gsub("_"," ",envVariableName)
						substr(envVariableName,1,1) = toupper(substr(envVariableName,1,1))				
						mtext(envVariableName, line=-1.8, cex=0.8, col="gray30")
						plot(envVariablesPlot[[i]], legend.only=T, add=T, col=colourScale, legend.width=0.5, alpha=1, legend.shrink=0.3, smallplot=c(0.35,0.87,0.14,0.17),
							 horizontal=T, legend.args=list(text="", cex=0.8, line=0.5, col.axis="gray30", col.lab="gray30", col="gray30"),
							 axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.2, tck=-0.8, line=0, mgp=c(0,0.4,0), col.tick="gray30", col.axis="gray30", col.lab="gray30", col="gray30"))
					}
				dev.off()
			}
	}

	# 1.4. Organising the environmental variables by sets

data_sources = read.csv("Environmental_data/Data_sources.csv", head=T, sep=";")
variables_sets = unique(data_sources[,"kind_of_variable"])
envVariables = list(); envVariableNames = list()
for (i in 1:length(variables_sets))
	{
		buffer1 = list()
		buffer2 = data_sources[,"variable_name_2"]
		buffer2 = buffer2[which(data_sources[,"kind_of_variable"]==variables_sets[i])]
		buffer2 = buffer2[which(buffer2!="")]
		buffer2 = buffer2[!buffer2%in%c("vaccination_in_china")]
		for (j in 1:length(buffer2))
			{		
				envVariableName = buffer2[j]; substr(envVariableName,1,1) = toupper(substr(envVariableName,1,1))				
				buffer1[[j]] = raster(paste0("Environmental_data/Prepared_rasters/",envVariableName,".tif"))
				if (grepl("population_density",envVariableName))
					{
						buffer1[[j]][] = log(buffer1[[j]][]+1)
						buffer2[j] = paste0(buffer2[j],"_(log)")
					}
			}
		envVariables[[i]] = buffer1; envVariableNames[[i]] = buffer2
	}
variables_sets = c(variables_sets, "set 5: new set for the wild birds"); buffer1 = list()
buffer2 = c("duck_population_density_2010",
			"extensive_chicken_population_density_2015",
			"intensive_chicken_population_density_2015",
			"human_population_density_2020",
			"urban_and_built_up_areas",
			"open_water_areas",
			"mixed_and_other_trees",
			"herbaceous_vegetation",
			"day_LST_amplitude_bi_annual",
			"NDVI_annual_mean",
			"EVI_annual_mean",
			"precipitation_February",
			"humidity_August")
for (i in 1:length(buffer2))
	{		
		envVariableName = buffer2[i]; substr(envVariableName,1,1) = toupper(substr(envVariableName,1,1))				
		buffer1[[length(buffer1)+1]] = raster(paste0("Environmental_data/Prepared_rasters/",envVariableName,".tif"))
		if (grepl("population_density",buffer2[i]))
			{
				buffer1[[i]][] = log(buffer1[[i]][]+1)
				buffer2[i] = paste0(buffer2[i],"_(log)")
			}
	}
envVariables[[length(envVariables)+1]] = buffer1; envVariableNames[[length(envVariableNames)+1]] = buffer2
variables_sets = c(variables_sets, "set 6: new set for the domestic birds"); buffer1 = list()
buffer2 = c("duck_population_density_2010",
			"extensive_chicken_population_density_2015",
			"intensive_chicken_population_density_2015",
			"human_population_density_2020",
			"urban_and_built_up_areas",
			"cultivated_and_managed_vegetation",
			"evergreen_broadleaf_trees",
			"deciduous_broadleaf_trees",
			"open_water_areas",
			"day_LST_annual_mean",
			"EVI_variance_tri_annual",
			"humidity_April")
for (i in 1:length(buffer2))
	{		
		envVariableName = buffer2[i]; substr(envVariableName,1,1) = toupper(substr(envVariableName,1,1))				
		buffer1[[length(buffer1)+1]] = raster(paste0("Environmental_data/Prepared_rasters/",envVariableName,".tif"))
		if (grepl("population_density",buffer2[i]))
			{
				buffer1[[i]][] = log(buffer1[[i]][]+1)
				buffer2[i] = paste0(buffer2[i],"_(log)")
			}
	}
envVariables[[length(envVariables)+1]] = buffer1; envVariableNames[[length(envVariableNames)+1]] = buffer2

# 2. Loading and inspecting the occurrence records

tab = read.csv("Occurrence_records/Dataset_06072023.csv", head=T, sep=";")
v_mask = raster::extract(mask, SpatialPoints(tab[,c("longitude","latitude")]))
tab = tab[!is.na(v_mask),]; serotype_summary = table(tab[,"serotype"]); serotypes = c("H5N1","H5")
dates = ymd(tab[,"report_date"]); decimal_dates = decimal_date(dates); cols = list()
cols[[1]] = rgb(250,165,33,187,maxColorValue=255); cols[[2]] = rgb(70,118,187,207,maxColorValue=255)
if (savingPlots)
	{
		pdf("H5_epi_curves_1_NEW.pdf", width=8, height=3); par(mfrow=c(2,1), oma=c(0,0,0.5,0.5), mar=c(1.3,1.3,0,0), lwd=0.2, col="gray30")
		hist(decimal_dates[which(grepl("H5N1",tab[,"serotype"]))], breaks=length(unique(dates))/7, xlim=c(2015,2023), ylim=c(0,680), ann=F, axes=F, col="gray80", border=NA)
		hist(decimal_dates[which(grepl("H5N1",tab[,"serotype"])&grepl("Wild",tab[,"bird_type"]))],
			 breaks=length(unique(dates))/7, xlim=c(2015,2023), ylim=c(0,680), add=T, ann=F, axes=F, col="gray65", border=NA)
		axis(side=1, pos=-30, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.030, col.axis="gray30", mgp=c(0,-0.1,0), at=c(2015:2024))
		axis(side=2, pos=2014.94, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.03, col.axis="gray30", mgp=c(0,0.18,0), at=c(0,200,400,600,800))
		title(ylab="H5N1 cases", cex.lab=0.7, mgp=c(0,0,0), col.lab="gray30")
		hist(decimal_dates[which(grepl("H5",tab[,"serotype"]))], breaks=length(unique(dates))/7, xlim=c(2015,2023), ylim=c(0,680), ann=F, axes=F, col="gray80", border=NA)
		hist(decimal_dates[which(grepl("H5",tab[,"serotype"])&grepl("Wild",tab[,"bird_type"]))], 
			 breaks=length(unique(dates))/7, xlim=c(2015,2023), ylim=c(0,680), add=T, ann=F, axes=F, col="gray65", border=NA)
		axis(side=1, pos=-30, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.030, col.axis="gray30", mgp=c(0,-0.1,0), at=c(2015:2024))
		axis(side=2, pos=2014.94, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.03, col.axis="gray30", mgp=c(0,0.18,0), at=c(0,200,400,600,800))
		title(ylab="H5Nx cases", cex.lab=0.7, mgp=c(0,0,0), col.lab="gray30")
		dev.off()
		pdf("H5_epi_curves_2_NEW.pdf", width=8, height=3); par(mfrow=c(2,1), oma=c(0,0,0.5,0.5), mar=c(1.5,1.3,0,0), lwd=0.2, col="gray30")
		hist(decimal_dates[which(grepl("H5",tab[,"serotype"])&grepl("Wild",tab[,"bird_type"]))],
			 breaks=length(unique(dates))/7, xlim=c(2015,2023), ylim=c(0,490), ann=F, axes=F, col=cols[[2]], border=NA)
		hist(decimal_dates[which(grepl("H5N1",tab[,"serotype"])&grepl("Wild",tab[,"bird_type"]))],
			 breaks=length(unique(dates))/7, xlim=c(2015,2023), ylim=c(0,490), add=T, ann=F, axes=F, col="white", border=NA)
		hist(decimal_dates[which(grepl("H5N1",tab[,"serotype"])&grepl("Wild",tab[,"bird_type"]))],
			 breaks=length(unique(dates))/7, xlim=c(2015,2023), ylim=c(0,490), add=T, ann=F, axes=F, col=cols[[1]], border=NA)
		axis(side=1, pos=-30, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.030, col.axis="gray30", mgp=c(0,-0.1,0), at=c(2015:2024))
		axis(side=2, pos=2014.94, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.03, col.axis="gray30", mgp=c(0,0.18,0), at=c(0,200,400,600))
		title(ylab="Wild bird cases", cex.lab=0.7, mgp=c(0,0,0), col.lab="gray30")
		hist(decimal_dates[which(grepl("H5",tab[,"serotype"])&grepl("Domestic",tab[,"bird_type"]))],
			 breaks=length(unique(dates))/7, xlim=c(2015,2023), ylim=c(0,490), ann=F, axes=F, col=cols[[2]], border=NA)
		hist(decimal_dates[which(grepl("H5N1",tab[,"serotype"])&grepl("Domestic",tab[,"bird_type"]))],
			 breaks=length(unique(dates))/7, xlim=c(2015,2023), ylim=c(0,490), add=T, ann=F, axes=F, col="white", border=NA)
		hist(decimal_dates[which(grepl("H5N1",tab[,"serotype"])&grepl("Domestic",tab[,"bird_type"]))],
			 breaks=length(unique(dates))/7, xlim=c(2015,2023), ylim=c(0,490), add=T, ann=F, axes=F, col=cols[[1]], border=NA)
		axis(side=1, pos=-30, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.030, col.axis="gray30", mgp=c(0,-0.1,0), at=c(2015:2024))
		axis(side=2, pos=2014.94, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.03, col.axis="gray30", mgp=c(0,0.18,0), at=c(0,200,400,600))
		title(ylab="Domestic bird cases", cex.lab=0.7, mgp=c(0,0,0), col.lab="gray30")
		dev.off()
	}
admin_pol = matrix(nrow=dim(tab)[1], ncol=1); colnames(admin_pol) = "admin_pol"
admins1_CHN = shapefile("Environmental_data/GADM_shapefiles/CHN_admin1.shp")
admins1_RUS = shapefile("Environmental_data/GADM_shapefiles/RUS_admin1.shp")
admin_pol[,"admin_pol"] = tab[,"country"]
for (i in 1:dim(tab)[1])
	{
		if (admin_pol[i,"admin_pol"] == "China")
			{
				polID = c()
				for (j in 1:length(admins1_CHN@polygons))
					{
						for (k in 1:length(admins1_CHN@polygons[[j]]@Polygons))
							{
						if (point.in.polygon(tab[i,"longitude"],tab[i,"latitude"],admins1_CHN@polygons[[j]]@Polygons[[k]]@coords[,1],admins1_CHN@polygons[[j]]@Polygons[[k]]@coords[,2]) == 1)
							{
								polID = c(polID,j)
							}}
					}
				if (length(polID) != 1) print(i)
				admin_pol[i,"admin_pol"] = admins1_CHN@data[polID,"NAME_1"]
			}
		if (admin_pol[i,"admin_pol"] == "Russian Federation")
			{
				polID = c()
				for (j in 1:length(admins1_RUS@polygons))
					{
						for (k in 1:length(admins1_RUS@polygons[[j]]@Polygons))
							{
						if (point.in.polygon(tab[i,"longitude"],tab[i,"latitude"],admins1_RUS@polygons[[j]]@Polygons[[k]]@coords[,1],admins1_RUS@polygons[[j]]@Polygons[[k]]@coords[,2]) == 1)
							{
								polID = c(polID,j)
							}}
					}
				if (length(polID) != 1) print(i)
				admin_pol[i,"admin_pol"] = admins1_RUS@data[polID,"NAME_1"]
			}
	}
tab = cbind(tab, admin_pol)
datasets1 = list(); datasets2 = list() # in the latest: only admins with at least 5 observations
dataset_IDs = c() # H5N1 and H5Nx, before and after the start of 2020, wild and domestic cases
for (i in 1:length(serotypes))
	{
		sub1 = tab[which(grepl(serotypes[i],tab[,"serotype"])),]; serotype = serotypes[i]
		dates_sub1 = ymd(sub1[,"report_date"]); decimal_dates_sub1 = decimal_date(dates_sub1)
		for (j in 1:2)
			{
				if (j == 1) { sub2 = sub1[which(decimal_dates_sub1<2020),]; period = "t1" }
				if (j == 2) { sub2 = sub1[which(decimal_dates_sub1>=2020),]; period = "t2" }
				for (k in 1:2)
					{
						if (k == 1) { sub3 = sub2[which(sub2[,"bird_type"]=="Wild"),]; host = "wild" }
						if (k == 2) { sub3 = sub2[which(sub2[,"bird_type"]=="Domestic"),]; host = "domestic" }
						datasets1[[length(datasets1)+1]] = sub3; sub4 = sub3[1,]
						lon_lat = paste0(sub3[1,"longitude"],"_",sub3[1,"latitude"])
						for (l in 2:dim(sub3)[1])
							{
								buffer = paste0(sub3[l,"longitude"],"_",sub3[l,"latitude"])
								if (!buffer%in%lon_lat)
									{
										sub4 = rbind(sub4, sub3[l,]); lon_lat = c(lon_lat, buffer)
									}
							}
						different_admins = table(sub4[,"admin_pol"])
						admins_to_keep = names(different_admins)[as.numeric(which(different_admins>=5))]						
						sub5 = sub3[which(sub3[,"admin_pol"]%in%admins_to_keep),]; datasets2[[length(datasets2)+1]] = sub5
						dataset_IDs = c(dataset_IDs, paste(serotype,period,host,sep="_"))
					}
			}
	}

# 3. Preparing the BRT dataframe for each dataset

if (!file.exists("BRT_dataframes.rds"))
	{
		dataframes = list()
		pseudoAbsences = list(); mask = raster("Environmental_data/Prepared_rasters/Mask_raster_08.tif")
		human_pop_density = raster("Environmental_data/Prepared_rasters/Human_population_density_2020.tif")
		human_pop_density_log = human_pop_density; human_pop_density_log[] = 1+log(human_pop_density_log[]+1)
		if (!file.exists("Environmental_data/GADM_combined/GADM_combined.shp"))
			{
				iso3s = read.csv("Country_ISO2_3.csv", head=T, sep=";")[,c("country","ISO3")]
				different_countries = unique(tab[,"country"]); different_countries = different_countries[order(different_countries)]
				for (i in 1:length(different_countries))
					{
						iso3 = iso3s[which(iso3s[,"country"]==different_countries[i]),"ISO3"]
						if (length(iso3) == 0) print(c(i,different_countries[i],iso3))
						if ((iso3 != "CHN")&(iso3 != "RUS"))
							{
								if (i == 1)
									{
										gadms = raster::getData("GADM", country=iso3, level=0)
										gadms@data[,"NAME_0"] = different_countries[i]
									}	else	{
										gadm = raster::getData("GADM", country=iso3, level=0)
										gadm@data[,"NAME_0"] = different_countries[i]
										gadms = rbind(gadms, gadm, makeUniqueIDs=T)
									}
							}
						if (iso3 == "CHN")
							{
								for (j in 1:length(admins1_CHN@polygons))
									{
										name1 = admins1_CHN@data[j,"NAME_1"]
										admin1 = subset(admins1_CHN, NAME_1==name1)
										admin1@data = admin1@data[,c("GID_1","NAME_1")]
										names(admin1) = c("GID_0","NAME_0")
										gadms = rbind(gadms, admin1, makeUniqueIDs=T)
									}
							}
						if (iso3 == "RUS")
							{
								for (j in 1:length(admins1_RUS@polygons))
									{
										name1 = admins1_RUS@data[j,"NAME_1"]
										admin1 = subset(admins1_RUS, NAME_1==name1)
										admin1@data = admin1@data[,c("GID_1","NAME_1")]
										names(admin1) = c("GID_0","NAME_0")
										gadms = rbind(gadms, admin1, makeUniqueIDs=T)
									}
							}
					}
				writeOGR(gadms, dsn="Environmental_data/GADM_combined/", layer="GADM_combined", driver="ESRI Shapefile")
			}	else	{
				gadms = shapefile("Environmental_data/GADM_combined/GADM_combined.shp")
			}
		for (i in 1:length(datasets2))
			{
				dataframe = cbind(rep(1,dim(datasets2[[i]])[1]), datasets2[[i]][,"longitude"], datasets2[[i]][,"latitude"]); buffer = c()
				admins_to_select = unique(datasets2[[i]][,"admin_pol"]); admins_to_select = admins_to_select[order(admins_to_select)]
				for (j in 1:length(admins_to_select))				
					{
						admin = subset(gadms, NAME_0%in%admins_to_select[j])
						background = mask(crop(human_pop_density_log, admin), admin)
						nberOfPresencePoints = length(which(datasets2[[i]][,"admin_pol"]==admins_to_select[j]))
						nonNAcells = length(which(!is.na(background[]))); N = 10*nberOfPresencePoints
						if (N > nonNAcells) N = nonNAcells
						PAs = seegSDM::bgSample(background, n=N, prob=T, replace=F, spatial=F)
							# sampling probability proportional to the log-transformed number of people living in each raster cell
						minimumDistances = rep(NA, dim(PAs)[1])
						for (k in 1:length(minimumDistances))
							{
								minimumDistances[k] = min(rdist.earth(cbind(PAs[k,1],PAs[k,2]), cbind(datasets2[[i]][,"longitude"],datasets2[[i]][,"latitude"]), miles=F))
							}
						PAs = PAs[which((minimumDistances>10)&(minimumDistances<1000)),] # discarding pseudo-absence points <10 or >1000 km of presence points are discarded
						if ((3*nberOfPresencePoints) < dim(PAs)[1])
							{
								PAs = PAs[sample(1:dim(PAs)[1],3*nberOfPresencePoints),]
							}
						buffer = rbind(buffer, PAs)
						if (showingPlots)
							{
								plot(background, box=F, axes=F, legend=F); points(PAs)
								points(datasets2[[i]][which(datasets2[[i]][,"admin_pol"]==admins_to_select[j]),c("longitude","latitude")], col="red")
							}
					}
				if (showingPlots)
					{
						plot(human_pop_density_log, box=F, axes=F, legend=F); points(buffer)
						points(datasets2[[i]][,c("longitude","latitude")], col="red")
					}
				dataframe = rbind(dataframe, cbind(rep(0,dim(buffer)[1]),buffer)); colnames(dataframe) = c("response","longitude","latitude")
				data_to_discard = which(is.na(raster::extract(mask,dataframe[,c("longitude","latitude")])))
				data_to_discard = data_to_discard[order(data_to_discard)]
				dataframe = dataframe[which(!c(1:dim(dataframe)[1])%in%data_to_discard),]		
				cellIDs = cellFromXY(mask, dataframe[,c("longitude","latitude")]); buffer = c()
				for (j in 1:length(unique(cellIDs))) # to keep only one presence or pseudo-absence point per raster cell 
					{								 # (priority = presence points):
						if (sum(cellIDs==unique(cellIDs)[j]) > 1)
							{
								tmp = dataframe[which(cellIDs==unique(cellIDs)[j]),]
								if (sum(tmp[,"response"]==1) == 0)
									{
										buffer = rbind(buffer, tmp[sample(1:dim(tmp)[1],1),])
									}
								if (sum(tmp[,"response"]==1) == 1)
									{
										buffer = rbind(buffer, tmp[which(tmp[,"response"]==1),])
									}
								if (sum(tmp[,"response"]==1) >= 2)
									{
										indices = which(tmp[,"response"]==1)
										buffer = rbind(buffer, tmp[sample(indices,1),])
									}
							}	else	{
								buffer = rbind(buffer, dataframe[which(cellIDs==unique(cellIDs)[j]),])
							}
					}
				dataframe = buffer
				buffer = matrix(nrow=dim(dataframe)[1], ncol=length(allEnvVariableNames))
				colnames(buffer) = allEnvVariableNames
				for (j in 1:dim(buffer)[2])
					{
						buffer[,j] = raster::extract(allEnvVariables[[j]], dataframe[,c("longitude","latitude")])
					}
				dataframes[[i]] = cbind(dataframe, buffer)
			}
		saveRDS(dataframes, "BRT_dataframes.rds")
	}
dataframes = readRDS("BRT_dataframes.rds")

# 4. All boosted regression trees (BRT) analyses

samplingPtsMinDist = function(observations, minDist=500, nberOfPoints=5)
	{
		indices = rep(NA, nberOfPoints)
		selection_list = list(1:nrow(observations)) 
  		indices[1] = sample(1:dim(observations)[1], 1)
		dists = list(spDistsN1(as.matrix(observations), as.matrix(observations[indices[1],]), longlat=T))
		for (i in 2:nberOfPoints)
			{
    			selection = which(dists[[(i-1)]] > minDist)
    			if (length(selection) == 0)
    				{
    					stop("Restarts the function with a smaller minimum distance")
					}
    			selection_list[[i]] = selection
    			test = table(unlist(selection_list))
    			indices_minDist = as.numeric(names(which(test==i)))
    			indices[i] = sample(indices_minDist, 1)   
				dists[[i]] = spDistsN1(as.matrix(observations), as.matrix(observations[indices[i],]), longlat=T)
			}
		return(indices)
	}
foldSelection = function(observations, selectedPoints)
	{
		fold_selection = sapply(1:nrow(observations), function(i) which.min(spDistsN1(as.matrix(selectedPoints), as.matrix(observations[i,]), longlat=T)))
		return(fold_selection)
	}
plottingCorrelogram = FALSE
if (plottingCorrelogram == TRUE)
	{
		for (i in 1:length(dataframes)) # distances from which the correlation values reach 0: ~2000, ~1750, 2000, 2000
			{
				buffer = dataframes[[i]][sample(1:dim(dataframes[[i]])[1],5000,replace=F),]
				correlogram = ncf::correlog(buffer[,"longitude"], buffer[,"latitude"], buffer[,"response"], na.rm=T, increment=100, resamp=0, latlon=T)
				pdf(paste0("Correlogram_",i,".pdf"), width=4.5, height=3); par(mar=c(2.2,2.2,1.5,1.5))
				plot(correlogram$mean.of.class[-1], correlogram$correlation[-1], ann=F, axes=F, lwd=0.2, cex=0.5, col=NA, ylim=c(-0.4,1.0), xlim=c(0,8500))
				abline(h=0, lwd=0.5, col="red", lty=2)
				points(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, cex=0.35, col="gray30")
				lines(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, col="gray30")
				axis(side=1, pos=-0.4, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,-0.05,0), at=seq(0,9000,1000))
				axis(side=2, pos=0, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.18,0), at=seq(-0.4,1,0.2))
				title(xlab="distance (km2)", cex.lab=0.7, mgp=c(0.3,0,0), col.lab="gray30")
				title(ylab="correlation", cex.lab=0.7, mgp=c(0.4,0,0), col.lab="gray30"); dev.off()
			}
	}
newAnalyses = TRUE; newAnalyses = FALSE; nberOfReplicates = 10 # one replicate = one folds partition
if (newAnalyses == TRUE) { for (t in 1:length(variables_sets)) { for (i in 2:length(datasets2)) {
		if (!file.exists(paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_CCV_SCV_AUCs.csv")))
			{
				AUCs = matrix(nrow=nberOfReplicates, ncol=3); colnames(AUCs) = c("CCV_AUC","SCV1_AUC","SCV2_AUC")
			}	else	{
				AUCs = read.csv(paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_CCV_SCV_AUCs.csv"), head=T)
			}
		mask = raster("Environmental_data/Prepared_rasters/Mask_raster_08.tif")
		data = dataframes[[i]]; envVariableColNames = colnames(data)[4:dim(data)[2]]
		variables_to_keep = envVariableColNames[which(envVariableColNames%in%envVariableNames[[t]])]
		data = data.frame(data[,c(1:3,which(colnames(data)%in%variables_to_keep))])
		colnames(data) = gsub("\\.","",colnames(data))
		variables_to_keep = gsub("\\(","",gsub("\\)","",variables_to_keep))
		theRanges = c(2000,2000)*1000 # distance in meters
		gbm.x = variables_to_keep
		gbm.y = "response"
		offset = NULL
		tree.complexity = 4 # "tc" = number of nodes in the trees
		learning.rate = 0.01 # "lr" = contribution of each tree to the growing model
		bag.fraction = 0.80 # proportion of data used to train a given tree
		site.weights = rep(1, dim(data)[1])
		var.monotone = rep(0, length(gbm.x))
		n.folds = 5
		prev.stratify = TRUE
		family = "bernoulli"
		n.trees = 100 # initial number of trees
		step.size = 10 # interval at which the predictive deviance is computed and logged
					  # (at each interval, the folds are successively used as test data set
					  # nd the remaining folds as training data sets to compute the deviance)
		max.trees = 10000 # maximum number of trees that will be considered
		tolerance.method = "auto"
		tolerance = 0.001
		plot.main = TRUE
		plot.folds = FALSE
		verbose = TRUE
		silent = FALSE
		keep.fold.models = FALSE
		keep.fold.vector = FALSE
		keep.fold.fit = FALSE
		showingFoldsPlot = FALSE
		brt_model_ccvs = list() # classic cross-validations (CCVs)
		brt_model_scv1 = list() # spatial cross-validations 1 (SCV1)
		brt_model_scv2 = list() # spatial cross-validations 2 (SCV2)
		for (j in 1:nberOfReplicates)
			{
				if (!file.exists(paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_CCV_replicate_",j,".rds")))
					{
						# BRT with classic (standard) cross-validation (CCV):
						pdf(file=paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_CCV_replicate_",j,".pdf"))
						n.folds = 5; n.trees = 100; fold.vector = NULL
						brt_model_ccvs[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
							var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
							verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
						dev.off()
						saveRDS(brt_model_ccvs[[j]], paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_CCV_replicate_",j,".rds"))
					}	else	{
						brt_model_ccvs[[j]] = readRDS(paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_CCV_replicate_",j,".rds"))
					}
				AUCs[j,1] = brt_model_ccvs[[j]]$cv.statistics$discrimination.mean # Mean test AUC (from the AUCs computed on each fold tested as test data in the CCV)
				if (!file.exists(paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_SCV1_replicate_",j,".rds")))
					{
						# BRT with spatial (geographic) cross-validation (SCV) based on the folds generation of Dhingra, Artois et al. (2016, eLife):
; 						worked = FALSE
						while (worked == FALSE)
							{
								trycatch = tryCatch(
									{
										n.folds = 5; n.trees = 100; c = 0; folds_with_similar_sizes = FALSE
										while (folds_with_similar_sizes == FALSE) # while loop to select a partition where the x folds gather at least
											{									  # a proportion = (1/(x+1)) of the total number of presence points
												data_presence = data[which(data[,"response"]==1),]; c = c+1; # print(c)
												points = samplingPtsMinDist(data_presence[,c("longitude","latitude")], minDist=500, nberOfPoints=n.folds)
												fold.vector = foldSelection(data[,c("longitude","latitude")], selectedPoints=data_presence[points,c("longitude","latitude")])
												fold.vector_presences = fold.vector[which(data[,"response"]==1)]
												counts = hist(fold.vector_presences, plot=F)$counts
												props = counts[which(counts > 0)]/sum(counts); print(round(props,2))
												if ((length(props) == n.folds)&(min(props) > (1/(2*n.folds)))) folds_with_similar_sizes = TRUE
											}
										if (showingFoldsPlot == TRUE)
											{
												par(mar=c(0,0,0,0), oma=c(0.0,3.6,0.0,0.0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
												cols = c("olivedrab3","tan3","steelblue3","orange1","tomato2","mediumseagreen")[fold.vector]
												plot(mask, col="gray90", useRaster=T, colNA=NA, box=F, axes=F, legend=F)
												pchs = c(16,3)[data[,"response"]+1]; cexs = c(0.25,0.5)[data[,"response"]+1]
												points(data[,c("longitude","latitude")], col=cols, pch=pchs, cex=cexs, lwd=0.7)
											}
										pdf(file=paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_SCV1_replicate_",j,".pdf"))
										brt_model_scv1[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
											var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
											verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit)
										dev.off()
										saveRDS(brt_model_scv1[[j]], paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_SCV1_replicate_",j,".rds"))
									},	error = function(cond) {
									},	finally = {
									})
								if (length(brt_model_scv1) == j) worked = TRUE
							}
					}	else	{
						brt_model_scv1[[j]] = readRDS(paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_SCV1_replicate_",j,".rds"))
					}
				AUCs[j,"SCV1_AUC"] = brt_model_scv1[[j]]$cv.statistics$discrimination.mean
					# Mean test AUC (from the AUCs computed on each fold tested as test data in the SCV)		
				if (!file.exists(paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_SCV2_replicate_",j,".rds")))
					{
						# BRT with spatial (geographic) cross-validation (SCV) based on the blocks generation of Valavi et al. (2019, MEE):
 						worked1 = FALSE
						while (worked1 == FALSE)
							{
								trycatch = tryCatch(
									{
										spdf = SpatialPointsDataFrame(data[c("longitude","latitude")], data[,c(1,4:dim(data)[2])], proj4string=crs(mask))
										n.folds = 5; n.trees = 100; worked2 = FALSE
										while (worked2 == FALSE)
											{
												trycatch = tryCatch(
													{
														myblocks = NULL
														# myblocks = spatialBlock(spdf, species="response", rasterLayer=mask, k=n.folds, theRange=theRanges[1], selection="random")
														myblocks = cv_spatial(x=spdf, column="response", r=mask, size=theRanges[1], k=n.folds, hexagon=F, selection="random")
													},	error = function(cond) {
													},	finally = {
													})
												if (!is.null(myblocks)) worked2 = TRUE
											}
										fold.vector = myblocks$folds_ids
										pdf(file=paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_SCV2_replicate_",j,".pdf"))
										brt_model_scv2[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
											var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
											verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
										dev.off()
										saveRDS(brt_model_scv2[[j]], paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_SCV2_replicate_",j,".rds"))
									},	error = function(cond) {
									},	finally = {
									})
								if (length(brt_model_scv2) == j) worked1 = TRUE
							}
					}	else	{
						brt_model_scv2[[j]] = readRDS(paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_SCV2_replicate_",j,".rds"))
					}
				AUCs[j,"SCV2_AUC"] = brt_model_scv2[[j]]$cv.statistics$discrimination.mean
					# Mean test AUC (from the AUCs computed on each fold tested as test data in the SCV)
			}
		write.csv(AUCs, paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_CCV_SCV_AUCs.csv"), row.names=F, quote=F)
	}}}
if (newAnalyses == TRUE) { for (i in 2:length(datasets2)) {
		SSBs = matrix(nrow=nberOfReplicates, ncol=3); colnames(SSBs) = c("CCV_SSB","SCV1_SSB","SCV2_SSB")		
		t = 1; n.folds = 5; theRanges = c(2000,2000)*1000
		for (j in 1:nberOfReplicates)
			{
				brt_model = readRDS(paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_CCV_replicate_",j,".rds"))
				data = brt_model$gbm.call$dataframe[,1:3]; fold.vector = rep(NA, dim(data)[1])
				for (k in 1:dim(data)[1]) fold.vector[k] = sample(1:n.folds,1); vS = 0
				for (k in 1:n.folds)
					{
						p = data[which((data[,"response"]==1)&(fold.vector==k)),2:3]
						a = data[which((data[,"response"]==0)&(fold.vector==k)),2:3]
						reference = data[which((data[,"response"]==1)&(fold.vector!=k)),2:3]
						SSB = ssb(p, a, reference); vS = vS + (SSB[1,"p"]/SSB[1,"a"])
					}
				SSBs[j,"CCV_SSB"] = vS/n.folds
			}
		for (j in 1:nberOfReplicates)
			{
				brt_model = readRDS(paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_CCV_replicate_",j,".rds"))
				data = brt_model$gbm.call$dataframe[,1:3]; fold.vector = rep(NA, dim(data)[1])
				folds_with_similar_sizes = FALSE; c = 0; vS = 0
				while (folds_with_similar_sizes == FALSE)
					{
						data_presence = data[which(data[,"response"]==1),]; c = c+1
						fivePoints = samplingPtsMinDist(data_presence[,2:3], minDist=200, nberOfPoints=n.folds)
						fold.vector = foldSelection(data[,2:3], selectedPoints=data_presence[fivePoints,2:3])
						fold.vector_presences = fold.vector[which(data[,"response"]==1)]
						counts = hist(fold.vector_presences, plot=F)$counts
						props = counts[which(counts > 0)]/sum(counts); print(round(props,2))
						if (min(props) > (1/(n.folds*2))) folds_with_similar_sizes = TRUE
					}
				for (k in 1:n.folds)
					{
						p = data[which((data[,"response"]==1)&(fold.vector==k)),2:3]
						a = data[which((data[,"response"]==0)&(fold.vector==k)),2:3]
						reference = data[which((data[,"response"]==1)&(fold.vector!=k)),2:3]
						SSB = ssb(p, a, reference); vS = vS + (SSB[1,"p"]/SSB[1,"a"])
					}
				SSBs[j,"SCV1_SSB"] = vS/n.folds
			}
		for (j in 1:nberOfReplicates)
			{
				brt_model = readRDS(paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_CCV_replicate_",j,".rds"))
				data = brt_model$gbm.call$dataframe; fold.vector = rep(NA, dim(data)[1])
				worked = FALSE; folds_with_similar_sizes = FALSE; c = 0; vS = 0
				spdf = SpatialPointsDataFrame(data[c("longitude","latitude")], data[,c(1,4:dim(data)[2])], proj4string=crs(mask))
				while (worked == FALSE)
					{
						trycatch = tryCatch(
							{
								myblocks = NULL
								# myblocks = spatialBlock(spdf, species="response", rasterLayer=nullRasters[[t]], k=n.folds, theRange=theRanges[1], selection="random")
								myblocks = cv_spatial(x=spdf, column="response", r=mask, size=theRanges[1], k=n.folds, hexagon=F, selection="random")
							},	error = function(cond) {
							},	finally = {
							})
						if (!is.null(myblocks)) worked = TRUE
					}
				for (k in 1:n.folds)
					{
						fold.vector = myblocks$folds_ids
						p = data[which((data[,"response"]==1)&(fold.vector==k)),2:3]
						a = data[which((data[,"response"]==0)&(fold.vector==k)),2:3]
						reference = data[which((data[,"response"]==1)&(fold.vector!=k)),2:3]
						SSB = ssb(p, a, reference); vS = vS + (SSB[1,"p"]/SSB[1,"a"])
					}
				SSBs[j,"SCV2_SSB"] = vS/n.folds
			}
		write.csv(SSBs, paste0("All_the_BRT_models/Dataset_",i,"_CCV_SCV_all_the_SSB_estimates.csv"), row.names=F, quote=F)
	}}

# 5. Comparing all the AUC support and SSB values

	# SSB = Dp/Da (Hijmans 2012, Ecology), where:
		# Dp = mean distance between testing presence sites and nearest training-presence sites
		# Da = mean distance between testing absence sites and nearest training-presence sites
		# --> SSB = 1 suggests there is no spatial sorting bias
		# --> SSB = 0 suggests extreme spatial sorting bias

SSB_values = list()
for (i in 2:length(datasets2))
	{
		tab = read.csv(paste0("All_the_BRT_models/Dataset_",i,"_CCV_SCV_all_the_SSB_estimates.csv"), head=T)
		for (j in 1:dim(tab)[2]) SSB_values[[length(SSB_values)+1]] = tab[,j]
		SSB_values[[length(SSB_values)+1]] = NA
		if (i != length(datasets2)) SSB_values[[length(SSB_values)+1]] = NA
	}
if (savingPlots)
	{
		pdf("All_SSB_values_NEW.pdf", width=8, height=1.75); par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1.0,2.0,1.0,0.0), lwd=0.2, col="gray30")
		boxplot(SSB_values, ann=F, axes=F, col="gray80", border="gray30", lwd=0.4, outline=F)
		axis(side=2, pos=-0.5, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.03, col.axis="gray30", mgp=c(0,0.18,0))
		title(ylab="SSB", cex.lab=0.7, mgp=c(0.6,0,0), col.lab="gray30")
		dev.off()
	}
AUC_values_SC1 = list(); AUC_values_SC2 = list()
for (i in 2:length(datasets2))
	{
		for (t in 1:4) # length(variables_sets)
			{
				tab = read.csv(paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_CCV_SCV_AUCs.csv"), head=T)
				AUC_values_SC1[[length(AUC_values_SC1)+1]] = tab[,"SCV1_AUC"]
				AUC_values_SC2[[length(AUC_values_SC2)+1]] = tab[,"SCV2_AUC"]
			}
		if (i != length(datasets2))
			{
				AUC_values_SC1[[length(AUC_values_SC1)+1]] = NA; AUC_values_SC2[[length(AUC_values_SC2)+1]] = NA
			}
	}
if (savingPlots)
	{
		pdf("SCV1_AUC_vals_NEW.pdf", width=8, height=1.75); par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1.0,2.0,1.0,0.0), lwd=0.2, col="gray30")
		boxplot(AUC_values_SC1, ann=F, axes=F, col="gray80", border="gray30", lwd=0.4, outline=F)
		axis(side=2, pos=-0.5, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.03, col.axis="gray30", mgp=c(0,0.18,0), at=seq(0.5,1.0,0.1))
		title(ylab="AUC", cex.lab=0.7, mgp=c(0.6,0,0), col.lab="gray30")
		dev.off()
		pdf("SCV2_AUC_vals_NEW.pdf", width=8, height=1.75); par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1.0,2.0,1.0,0.0), lwd=0.2, col="gray30")
		boxplot(AUC_values_SC2, ann=F, axes=F, col="gray80", border="gray30", lwd=0.4, outline=F)
		axis(side=2, pos=-0.5, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.03, col.axis="gray30", mgp=c(0,0.18,0), at=seq(0.5,1.0,0.1))
		title(ylab="AUC", cex.lab=0.7, mgp=c(0.6,0,0), col.lab="gray30")
		dev.off()
	}

# 6. Computing the relative influence (RI) values

selectedModel = "SCV1"; selected_variables_sets = c(2,4,2,4,2,4,2,4)
dataset_names = c("H5N1_t1_wild","H5N1_t1_domestic","H5N1_t2_wild","H5N1_t2_domestic",
				  "H5Nx_t1_wild","H5Nx_t1_domestic","H5Nx_t2_wild","H5Nx_t2_domestic")
if (writingFiles)
	{
		for (t in 1:length(variables_sets))
			{
				RIs = matrix(nrow=length(envVariableNames[[t]]), ncol=length(datasets2))
				row.names(RIs) = envVariableNames[[t]]; colnames(RIs) = dataset_names
				for (i in 2:length(datasets2))
					{
						tmp = matrix(nrow=length(envVariableNames[[t]]), ncol=nberOfReplicates)
						for (j in 1:nberOfReplicates)
							{
								brt_model = readRDS(paste0("All_the_BRT_models/Variables_set_",t,"_dataset_",i,"_",selectedModel,"_replicate_",j,".rds"))
								for (k in 1:length(envVariableNames[[t]]))
									{
										tmp[k,j] = summary(brt_model)[gsub("\\(","",gsub("\\)","",envVariableNames[[t]][k])),"rel.inf"]
									}
							}
						for (j in 1:length(envVariableNames[[t]]))
							{
								RIs[j,i] = round(mean(tmp[j,]),1)
							}
					}
				write.csv(RIs, paste0("RI_values_set_",t,".csv"), quote=F)
			}
	}
if (writingFiles)
	{
		RIs = matrix(nrow=length(envVariableNames[[selected_variables_sets[3]]]), ncol=3)
		row.names(RIs) = envVariableNames[[selected_variables_sets[3]]]; ii = 0
		colnames(RIs) = c("H5N1_t2_wild","H5Nx_t1_wild","H5Nx_t2_wild")
		for (i in c(3,5,7))
			{
				tmp = matrix(nrow=length(envVariableNames[[selected_variables_sets[i]]]), ncol=nberOfReplicates); ii = ii+1
				for (j in 1:nberOfReplicates)
					{
						brt_model = readRDS(paste0("All_the_BRT_models/Variables_set_",selected_variables_sets[i],"_dataset_",i,"_",selectedModel,"_replicate_",j,".rds"))
						for (k in 1:length(envVariableNames[[selected_variables_sets[i]]]))
							{
								tmp[k,j] = summary(brt_model)[gsub("\\(","",gsub("\\)","",envVariableNames[[selected_variables_sets[i]]][k])),"rel.inf"]
							}
					}
				for (j in 1:length(envVariableNames[[selected_variables_sets[i]]]))
					{
						qts = round(quantile(tmp[j,], c(0.25,0.75)),1)
						RIs[j,ii] = paste0(round(median(tmp[j,]),1)," % [",qts[1],"-",qts[2],"]")
					}
			}
		write.csv(RIs, paste0("Table_S1_RIs_W.csv"), quote=F)
		RIs = matrix(nrow=length(envVariableNames[[selected_variables_sets[2]]]), ncol=4)
		row.names(RIs) = envVariableNames[[selected_variables_sets[2]]]; ii = 0
		colnames(RIs) = c("H5N1_t1_domestic","H5N1_t2_domestic","H5Nx_t1_domestic","H5Nx_t2_domestic")
		for (i in c(2,4,6,8))
			{
				tmp = matrix(nrow=length(envVariableNames[[selected_variables_sets[i]]]), ncol=nberOfReplicates); ii = ii+1
				for (j in 1:nberOfReplicates)
					{
						brt_model = readRDS(paste0("All_the_BRT_models/Variables_set_",selected_variables_sets[i],"_dataset_",i,"_",selectedModel,"_replicate_",j,".rds"))
						for (k in 1:length(envVariableNames[[selected_variables_sets[i]]]))
							{
								tmp[k,j] = summary(brt_model)[gsub("\\(","",gsub("\\)","",envVariableNames[[selected_variables_sets[i]]][k])),"rel.inf"]
							}
					}
				for (j in 1:length(envVariableNames[[selected_variables_sets[i]]]))
					{
						qts = round(quantile(tmp[j,], c(0.25,0.75)),1)
						RIs[j,ii] = paste0(round(median(tmp[j,]),1)," % [",qts[1],"-",qts[2],"]")
					}
			}
		write.csv(RIs, paste0("Table_S1_RIs_D.csv"), quote=F)
	}

# 7. Comparing the predictive capacities (with AUC)

newAnalyses = FALSE
if (newAnalyses == TRUE)
	{
		AUC_comparison = matrix(nrow=nberOfReplicates, ncol=9)
		colnames(AUC_comparison) = c("AUC_t1_t1_H5N1_domestic","AUC_t2_t2_H5N1_domestic","AUC_t1_t2_H5N1_domestic",
									 "AUC_t1_t1_H5Nx_wild","AUC_t2_t2_H5Nx_wild","AUC_t1_t2_H5Nx_wild",
									 "AUC_t1_t1_H5Nx_domestic","AUC_t2_t2_H5Nx_domestic","AUC_t1_t2_H5Nx_domestic")
		for (j in 1:nberOfReplicates)
			{
				c = 0
				for (i in c(2,5,6))
					{
						I = i; brt_model_1 = readRDS(paste0("All_the_BRT_models/Variables_set_",selected_variables_sets[I],"_dataset_",I,"_",selectedModel,"_replicate_",j,".rds"))
						I = i+2; brt_model_2 = readRDS(paste0("All_the_BRT_models/Variables_set_",selected_variables_sets[I],"_dataset_",I,"_",selectedModel,"_replicate_",j,".rds"))
						df = brt_model_1$gbm.call$dataframe; responses = df$response
						df = brt_model_1$gbm.call$dataframe; data = df[,4:dim(df)[2]]
						n.trees = brt_model_1$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(brt_model_1, data, n.trees, type, single.tree)
						c = c+1; AUC_comparison[j,c] = gbm.roc.area(obs=responses, pred=prediction)
						df = brt_model_2$gbm.call$dataframe; responses = df$response
						df = brt_model_2$gbm.call$dataframe; data = df[,4:dim(df)[2]]
						n.trees = brt_model_2$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(brt_model_2, data, n.trees, type, single.tree)
						c = c+1; AUC_comparison[j,c] = gbm.roc.area(obs=responses, pred=prediction)
						df = brt_model_2$gbm.call$dataframe; responses = df$response
						df = brt_model_2$gbm.call$dataframe; data = df[,4:dim(df)[2]]
						n.trees = brt_model_1$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(brt_model_1, data, n.trees, type, single.tree)
						c = c+1; AUC_comparison[j,c] = gbm.roc.area(obs=responses, pred=prediction)
					}
			}
		write.csv(round(AUC_comparison,3), "AUC_comparison_NEW.csv", row.names=F, quote=F)
		for (i in 1:dim(AUC_comparison)[2])
			{
				medianV = median(AUC_comparison[,i]); minV = min(AUC_comparison[,i]); maxV = max(AUC_comparison[,i])
				cat(paste0("\t\t\t\t# ",colnames(AUC_comparison)[i],": ",round(medianV,2)," [",round(minV,2),"-",round(maxV,2),"]\n"))				
			}
				# AUC_t1_t1_H5N1_domestic:	0.89 [0.89-0.90]
				# AUC_t2_t2_H5N1_domestic:	0.79 [0.77-0.80]
				# AUC_t1_t2_H5N1_domestic:	0.70 [0.70-0.71]
				# AUC_t1_t1_H5Nx_wild: 		0.84 [0.83-0.84]
				# AUC_t2_t2_H5Nx_wild:		0.82 [0.81-0.82]
				# AUC_t1_t2_H5Nx_wild:		0.77 [0.77-0.77]
				# AUC_t1_t1_H5Nx_domestic:	0.83 [0.82-0.83]
				# AUC_t2_t2_H5Nx_domestic:	0.78 [0.78-0.80]
				# AUC_t1_t2_H5Nx_domestic:	0.75 [0.74-0.76]
	}

# 8. Comparing the response curves and RI values

newAnalyses = FALSE
if (newAnalyses == TRUE)
	{
		envVariableValues_list = list()
		for (i in 2:length(datasets2))
			{
				data = dataframes[[i]]; data = data[which(data[,"response"]==1),]
				envVariableValues = matrix(nrow=3, ncol=length(envVariableNames[[selected_variables_sets[i]]]))
				row.names(envVariableValues) = c("median","minV","maxV")
				colnames(envVariableValues) = envVariableNames[[selected_variables_sets[i]]]
				for (j in 1:length(envVariableNames[[selected_variables_sets[i]]]))
					{
						minV = min(data[,envVariableNames[[selected_variables_sets[i]]][j]], na.rm=T)
						maxV = max(data[,envVariableNames[[selected_variables_sets[i]]][j]], na.rm=T)
						medianV = median(data[,envVariableNames[[selected_variables_sets[i]]][j]], na.rm=T)
						envVariableValues[,j] = cbind(medianV, minV, maxV)
					}
				envVariableValues_list[[i]] = envVariableValues
			}
		indices = list(); indices[[1]] = c(3,5,7); indices[[2]] = c(2,4,6,8)
		for (h in 1:length(indices))
			{
				predictions1 = list(); dfs1 = list()
				minMaxYValues = matrix(nrow=length(datasets2), ncol=2)
				minMaxXValues = matrix(nrow=13, ncol=2)
				for (i in indices[[h]])
					{
						predictions2 = list(); dfs2 = list()
						for (j in 1:length(envVariableNames[[selected_variables_sets[i]]]))
							{
								predictions3 = list(); valuesInterval = 0.1; valuesInterval = (envVariableValues_list[[i]]["maxV",j]-envVariableValues_list[[i]]["minV",j])/100
								df = data.frame(matrix(nrow=length(seq(envVariableValues_list[[i]]["minV",j],envVariableValues_list[[i]]["maxV",j],valuesInterval)),
													   ncol=length(envVariableNames[[selected_variables_sets[i]]]))); colnames(df) = envVariableNames[[selected_variables_sets[i]]]
								for (k in 1:length(envVariableNames[[selected_variables_sets[i]]]))
									{
										valuesInterval = 0.1; valuesInterval = (envVariableValues_list[[i]]["maxV",k]-envVariableValues_list[[i]]["minV",k])/100
										if (j == k) df[,envVariableNames[[selected_variables_sets[i]]][k]] = seq(envVariableValues_list[[i]]["minV",k],envVariableValues_list[[i]]["maxV",k],valuesInterval)
										if (j != k) df[,envVariableNames[[selected_variables_sets[i]]][k]] = rep(envVariableValues_list[[i]]["median",k],dim(df)[1])
									}
								colnames(df) = gsub("\\(","",gsub("\\)","",envVariableNames[[selected_variables_sets[i]]]))
								for (k in 1:nberOfReplicates)
									{
										brt_model = readRDS(paste0("All_the_BRT_models/Variables_set_",selected_variables_sets[i],"_dataset_",i,"_",selectedModel,"_replicate_",k,".rds"))
										n.trees = brt_model$gbm.call$best.trees; type = "response"; single.tree = FALSE
										prediction = predict.gbm(brt_model, newdata=df, n.trees, type, single.tree)
										if ((j == 1)&(k == 1))
											{
												minMaxYValues[i,1] = min(prediction); minMaxYValues[i,2] = max(prediction)
											}	else	{
												if (j <= 13)
													{
														if (minMaxYValues[i,1] > min(prediction)) minMaxYValues[i,1] = min(prediction)
														if (minMaxYValues[i,2] < max(prediction)) minMaxYValues[i,2] = max(prediction)
													}
											}
										if ((i == indices[[h]][1])&(k == 1))
											{
												minMaxXValues[j,1] = min(df[,gsub("\\(","",gsub("\\)","",envVariableNames[[selected_variables_sets[i]]][j]))])
												minMaxXValues[j,2] = max(df[,gsub("\\(","",gsub("\\)","",envVariableNames[[selected_variables_sets[i]]][j]))])
											}	else	{
												if (j <= 13)
													{
														if (minMaxXValues[j,1] > min(df[,gsub("\\(","",gsub("\\)","",envVariableNames[[selected_variables_sets[i]]][j]))]))
															{
																minMaxXValues[j,1] = min(df[,gsub("\\(","",gsub("\\)","",envVariableNames[[selected_variables_sets[i]]][j]))])
															}
														if (minMaxXValues[j,2] < max(df[,gsub("\\(","",gsub("\\)","",envVariableNames[[selected_variables_sets[i]]][j]))]))
															{
																minMaxXValues[j,2] = max(df[,gsub("\\(","",gsub("\\)","",envVariableNames[[selected_variables_sets[i]]][j]))])
															}
													}
											}
										predictions3[[k]] = prediction
									}
								predictions2[[j]] = predictions3; dfs2[[j]] = df
							}
						predictions1[[i]] = predictions2; dfs1[[i]] = dfs2
					}
				selected_variables_wild_birds = c("deciduous_broadleaf_trees",
												  "mixed_and_other_trees",
												  "herbaceous_vegetation",
												  "cultivated_and_managed_vegetation",
												  "urban_and_built_up_areas",
												  "open_water_areas",
												  "distance_to_water")
				selected_variables_domestic_birds = c("duck_population_density_2010_(log)",
												  "extensive_chicken_population_density_2015_(log)",
												  "intensive_chicken_population_density_2015_(log)",
												  "human_population_density_2020_(log)",
												  "cultivated_and_managed_vegetation",
												  "open_water_areas",
												  "day_LST_annual_mean")
				RIs = read.csv(paste0("RI_values_set_",selected_variables_sets[i],".csv"), head=T)
				pdf(paste0("Response_curves_NEW",h,".pdf"), width=7.5, height=3.5)
				par(mfrow=c(4,7), oma=c(0.5,0.5,1,1), mar=c(1.2,1.2,0.2,0.2), lwd=0.2, col="gray30")
				for (i in indices[[h]])
					{
						if (h == 1) selected_variables = selected_variables_wild_birds
						if (h == 2) selected_variables = selected_variables_domestic_birds
						for (j in 1:length(envVariableNames[[selected_variables_sets[i]]]))
							{
								if (envVariableNames[[selected_variables_sets[i]]][j]%in%selected_variables)
									{
								for (k in 1:length(predictions1[[i]][[j]]))
									{
										if (k == 1)
											{
												plot(dfs1[[i]][[j]][,gsub("\\(","",gsub("\\)","",envVariableNames[[selected_variables_sets[i]]][j]))], predictions1[[i]][[j]][[k]],
													 col="gray30", ann=F, axes=F, lwd=0.2, type="l", xlim=c(minMaxXValues[j,1],minMaxXValues[j,2]), 
													 ylim=c(minMaxYValues[i,1],minMaxYValues[i,2]))
											}	else	{
												lines(dfs1[[i]][[j]][,gsub("\\(","",gsub("\\)","",envVariableNames[[selected_variables_sets[i]]][j]))], predictions1[[i]][[j]][[k]], 
													  col="gray30", lwd=0.2)
											}
									}
								box(lwd=0.2, col="gray30"); mtext(paste0(RIs[j,i]," %"), lin=-2, cex=0.4, col="gray30")
								axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0, tck=-0.060, col.axis="gray30", mgp=c(0,0.00,0))
								axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0, tck=-0.060, col.axis="gray30", mgp=c(0,0.23,0))
								if (j == 1) title(ylab="predicted values", cex.lab=0.7, mgp=c(1.3,0,0), col.lab="gray30")
								if (i == indices[[h]][length(indices[[h]])])
									{
										title(xlab=gsub("_"," ",envVariableNames[[selected_variables_sets[i]]][j]), cex.lab=0.7, mgp=c(0.9,0,0), col.lab="gray30")				
									}
									}
							}
					}
				dev.off()
			}
	}

# 9. Projections of H5 ecological suitabilities

predictions1 = list()
for (i in 2:length(datasets2))
	{
		if (!file.exists(paste0("All_the_BRT_models/Variables_set_",selected_variables_sets[i],"_dataset_",i,"_",selectedModel,"_average.asc")))
			{
				buffer = envVariables[[selected_variables_sets[i]]]
				rast_NA = buffer[[1]]; rast_NA[!is.na(rast_NA)] = 0
				for (j in 2:length(buffer)) rast_NA[is.na(buffer[[j]][])] = NA
				for (j in 1:length(buffer)) buffer[[j]][is.na(rast_NA[])] = NA
				predictions2 = list(); rasters_stack = stack(buffer)
				df = as.data.frame(rasters_stack); not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
				colnames(newdata) = gsub("\\(","",gsub("\\)","",envVariableNames[[selected_variables_sets[i]]]))
				for (j in 1:nberOfReplicates)
					{
						brt_model = readRDS(paste0("All_the_BRT_models/Variables_set_",selected_variables_sets[i],"_dataset_",i,"_",selectedModel,"_replicate_",j,".rds"))		
						n.trees = brt_model$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(brt_model, newdata, n.trees, type, single.tree)
						rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction; predictions2[[j]] = rast
					}
				rasts = stack(predictions2); predictions1[[i]] = mean(rasts)
				writeRaster(predictions1[[i]], paste0("All_the_BRT_models/Variables_set_",selected_variables_sets[i],"_dataset_",i,"_",selectedModel,"_average.asc"))
			}
	}
for (i in 2:length(datasets2))
	{	
		predictions1[[i]] = crop(raster(paste0("All_the_BRT_models/Variables_set_",selected_variables_sets[i],"_dataset_",i,"_",selectedModel,"_average.asc")), extent(-170,180,-90,90))
	}
colourScale = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(14)[2:12])[c(1:4,7:8,8,9,9,10)]
mask_cropped = crop(mask, extent(-170,180,-90,90)); plottingPresencePoints = TRUE; plottingPresencePoints = FALSE
coast_lines = gSimplify(crop(shapefile("Environmental_data/NaturalEarth_files/Coastline_borders.shp"), extent(-170,180,-56,90)), 0.01)
dataset_names = c("H5N1\nwild birds\n< 2020","H5N1\ndomestic birds\n< 2020","H5N1\nwild birds\n> 2020","H5N1\ndomestic birds\n> 2020",
				  "H5Nx\nwild birds\n< 2020","H5Nx\ndomestic birds\n< 2020","H5Nx\nwild birds\n> 2020","H5Nx\ndomestic birds\n> 2020")
pdf(paste0("All_predictions_NEW.pdf"), width=7.5, height=6.0) # dev.new(width=7.5, height=6.0)
par(mfrow=c(4,2), oma=c(1,1,1,1), mar=c(0,0,0,0), mgp=c(0,0,0), lwd=0.2, col="gray30"); plot.new()
for (i in c(3,5,7,2,4,6,8))
	{
		if (plottingPresencePoints)
			{
				plot(mask_cropped, ann=F, axes=F, box=F, legend=F, col="gray90")
				plot(coast_lines, lwd=0.1, col="gray60", add=T)
				colourP = rgb(222,67,39,200,maxColorValue=255); colourA = rgb(77,77,77,30,maxColorValue=255)
				points(dataframes[[i]][which(dataframes[[i]][,"response"]==0),c("longitude","latitude")], pch=16, cex=0.2, col=colourA)
				points(dataframes[[i]][which(dataframes[[i]][,"response"]==1),c("longitude","latitude")], pch=16, cex=0.2, col=colourP)
				# points(datasets2[[i]][,c("longitude","latitude")], pch=16, cex=0.2, col=colourP)
			}	else	{
				if (i %in% c(1,3,5,7)) { minValue = 0.05; maxValue = 0.50 }
				if (i %in% c(2,4,6,8)) { minValue = 0.05; maxValue = 0.50 }
				prediction = predictions1[[i]]
				prediction[which(prediction[]<minValue)] = minValue
				prediction[which(prediction[]>maxValue)] = maxValue
				index1 = (((min(prediction[],na.rm=T)-minValue)/(maxValue-minValue))*10)+1
				index2 = (((max(prediction[],na.rm=T)-minValue)/(maxValue-minValue))*10)+1
				cols = colourScale[index1:index2]
				plot(prediction, ann=F, axes=F, box=F, legend=F, col=cols)
				plot(coast_lines, lwd=0.1, col="gray60", add=T)
			}
		mtext(dataset_names[i], at=-170, line=-10, adj=0, cex=0.6, col="gray30")
		if ((!plottingPresencePoints)&(i %in% c(8)))
			{
				legend_raster = raster(as.matrix(c(minValue,maxValue)))
				plot(legend_raster, legend.only=T, add=T, col=colourScale, legend.width=0.5, alpha=1, legend.shrink=0.3, smallplot=c(0.910,0.925,0.00,0.55),
					 horizontal=F, legend.args=list(text="", cex=0.8, line=0.5, col.axis="gray30", col.lab="gray30", col="gray30"),
					 axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.2, tck=-0.8, line=0, mgp=c(0,0.6,0), col.tick="gray30", col.axis="gray30", col.lab="gray30", col="gray30"))
			}
	}
dev.off()
pdf(paste0("Comparisons_NEW1.pdf"), width=7.5, height=6.0) # dev.new(width=7.5, height=6.0)
par(mfrow=c(4,2), oma=c(1,1,1,1), mar=c(0,0,0,0), mgp=c(0,0,0), lwd=0.2, col="gray30")
for (i in c(3,5,7))
	{
		plot(mask_cropped, ann=F, axes=F, box=F, legend=F, col="gray90")
		plot(coast_lines, lwd=0.1, col="gray60", add=T)
		colourP = rgb(222,67,39,200,maxColorValue=255); colourA = rgb(77,77,77,30,maxColorValue=255)
		points(dataframes[[i]][which(dataframes[[i]][,"response"]==0),c("longitude","latitude")], pch=16, cex=0.2, col=colourA)
		points(dataframes[[i]][which(dataframes[[i]][,"response"]==1),c("longitude","latitude")], pch=16, cex=0.2, col=colourP)
		mtext(dataset_names[i], at=-170, line=-10, adj=0, cex=0.6, col="gray30")
		if (i %in% c(1,3,5,7)) { minValue = 0.05; maxValue = 0.50 }
		if (i %in% c(2,4,6,8)) { minValue = 0.05; maxValue = 0.50 }
		prediction = predictions1[[i]]
		prediction[which(prediction[]<minValue)] = minValue
		prediction[which(prediction[]>maxValue)] = maxValue
		index1 = (((min(prediction[],na.rm=T)-minValue)/(maxValue-minValue))*10)+1
		index2 = (((max(prediction[],na.rm=T)-minValue)/(maxValue-minValue))*10)+1
		cols = colourScale[index1:index2]
		plot(prediction, ann=F, axes=F, box=F, legend=F, col=cols)
		plot(coast_lines, lwd=0.1, col="gray60", add=T)
		if (i == 7)
			{
				legend_raster = raster(as.matrix(c(minValue,maxValue)))
				plot(legend_raster, legend.only=T, add=T, col=colourScale, legend.width=0.5, alpha=1, legend.shrink=0.3, smallplot=c(0.910,0.925,0.00,0.55),
					 horizontal=F, legend.args=list(text="", cex=0.8, line=0.5, col.axis="gray30", col.lab="gray30", col="gray30"),
					 axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.2, tck=-0.8, line=0, mgp=c(0,0.6,0), col.tick="gray30", col.axis="gray30", col.lab="gray30", col="gray30"))
			}
	}
dev.off()
pdf(paste0("Comparisons_NEW2.pdf"), width=7.5, height=6.0) # dev.new(width=7.5, height=6.0)
par(mfrow=c(4,2), oma=c(1,1,1,1), mar=c(0,0,0,0), mgp=c(0,0,0), lwd=0.2, col="gray30")
for (i in c(2,4,6,8))
	{
		plot(mask_cropped, ann=F, axes=F, box=F, legend=F, col="gray90")
		plot(coast_lines, lwd=0.1, col="gray60", add=T)
		colourP = rgb(222,67,39,200,maxColorValue=255); colourA = rgb(77,77,77,30,maxColorValue=255)
		points(dataframes[[i]][which(dataframes[[i]][,"response"]==0),c("longitude","latitude")], pch=16, cex=0.2, col=colourA)
		points(dataframes[[i]][which(dataframes[[i]][,"response"]==1),c("longitude","latitude")], pch=16, cex=0.2, col=colourP)
		mtext(dataset_names[i], at=-170, line=-10, adj=0, cex=0.6, col="gray30")
		if (i %in% c(1,3,5,7)) { minValue = 0.05; maxValue = 0.50 }
		if (i %in% c(2,4,6,8)) { minValue = 0.05; maxValue = 0.50 }
		prediction = predictions1[[i]]
		prediction[which(prediction[]<minValue)] = minValue
		prediction[which(prediction[]>maxValue)] = maxValue
		index1 = (((min(prediction[],na.rm=T)-minValue)/(maxValue-minValue))*10)+1
		index2 = (((max(prediction[],na.rm=T)-minValue)/(maxValue-minValue))*10)+1
		cols = colourScale[index1:index2]
		plot(prediction, ann=F, axes=F, box=F, legend=F, col=cols)
		plot(coast_lines, lwd=0.1, col="gray60", add=T)
		if (i == 8)
			{
				legend_raster = raster(as.matrix(c(minValue,maxValue)))
				plot(legend_raster, legend.only=T, add=T, col=colourScale, legend.width=0.5, alpha=1, legend.shrink=0.3, smallplot=c(0.910,0.925,0.00,0.55),
					 horizontal=F, legend.args=list(text="", cex=0.8, line=0.5, col.axis="gray30", col.lab="gray30", col="gray30"),
					 axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.2, tck=-0.8, line=0, mgp=c(0,0.6,0), col.tick="gray30", col.axis="gray30", col.lab="gray30", col="gray30"))
			}
	}
dev.off()

