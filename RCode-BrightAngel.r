##### Bright Angel trophic cascade analysis #####
## Last updated 1 July 2019 by J.D. Muehlbauer


##### Set up workspace ##### 

## Load/install requisite packages
source('https://github.com/jmuehlbauer-usgs/R-packages/blob/master/packload.r?raw=TRUE')
packload(c('devtools', 'lubridate', 'plots', 'bugR'))

## Set working directory (depends on whether Jeff or Megan)
if(Sys.info()[6]=='jmuehlbauer'){
	setwd('C:/Users/jmuehlbauer/Documents/Projects/BrightAngel')
} else{
	setwd('C:/Users/mdaubert/Documents/Bright Angel/R')
}


##### Read in data #####

## Read and write drift and benthic data from GCMRC
	## Note: Only need to run this function once. Then only if the data need updating for some reason
BAreadwrite <- function(){
	## Get all raw data from database
	dbdat <- 'https://raw.githubusercontent.com/jmuehlbauer-usgs/Database/master/'
	files <- c('DriftSample', 'DriftSpecimen', 'tbl_BenthicSample', 
		'tbl_BenthicSpecimen', 'SpeciesList')
	gcmrc <- lapply(paste0(dbdat, files, '.csv'), read.csv)
		names(gcmrc) <- c('dsamp', 'dspec', 'bsamp', 'bspec', 'sppl')
	## Format dates and times for use
	gcmrc1 <- gcmrc
		gcmrc1$dsamp$Date <- as.Date(gcmrc$dsamp$Date, format = '%m/%d/%Y')
		gcmrc1$bsamp$Date <- as.Date(gcmrc$bsamp$SampleDate, format = '%m/%d/%Y')
		gcmrc1$bsamp$ProcessDate <- as.Date(gcmrc$bsamp$DateProcessed, format = '%m/%d/%Y')
	## Add barcodes to benthic specimens
	gcmrc1$bspec$BarcodeID <- gcmrc1$bsamp[match(gcmrc1$bspec$SampleID, gcmrc1$bsamp$SampleID), 'BarcodeID']
	## Subset to just the samples of interest (Bright Angel in 2016-January 2017)
	gcmrc2 <- gcmrc1
		gcmrc2$dsamp <- gcmrc1$dsamp[gcmrc1$dsamp$Reach == 'BrightAngel' & (year(gcmrc1$dsamp$Date) == 2016 | 
			(year(gcmrc1$dsamp$Date) == 2017 & month(gcmrc1$dsamp$Date) == 1)),]
			gcmrc2$dsamp <- droplevels(gcmrc2$dsamp)
		gcmrc2$bsamp <- gcmrc1$bsamp[gcmrc1$bsamp$River == 'Bright Angel' & (year(gcmrc1$bsamp$Date) == 2016 | 
			(year(gcmrc1$bsamp$Date) == 2017 & month(gcmrc1$bsamp$Date) == 1)),]
			gcmrc2$bsamp <- droplevels(gcmrc2$bsamp)
	## Format benthic sample data to parallel drift, rename or drop columns
	gcmrc2$bsamp$Region <- 'GrandCanyon'
	gcmrc2$bsamp$Reach <- 'BrightAngel'
	gcmrc2$bsamp$SampleNumber <- gcmrc2$bsamp$DatasheetSampleNo
	gcmrc2$bsamp$Depth <- gcmrc2$bsamp$SampleDepth
	gcmrc3 <- gcmrc2
	gcmrc3$bsamp <- gcmrc2$bsamp[,c('BarcodeID', 'TripID', 'Region', 'Reach', 'Date', 'SampleNumber', 
		'RiverMile', 'Depth', 'GearID', 'SampleArea', 'EntererSample', 'Processor', 
		'ProcessDate', 'ProcessTime', 'EntererSpecimen', 'Notes')]
		gcmrc3$bsamp$TripID <- ifelse(month(gcmrc3$bsamp$Date) == 6, 'BA20160608',
			ifelse(month(gcmrc3$bsamp$Date) == 11, 'BA20161108', 
			ifelse(month(gcmrc3$bsamp$Date) == 1, 'BA20170118', 'BA20160831')))
		gcmrc3$bsamp$TripID <- as.factor(gcmrc3$bsamp$TripID)
	## Format benthic specimen data to include barcodes and total counts, rename or drop columns
	gcmrc3$bspec$CountC <- rowSums(gcmrc3$bspec[,6:37], na.rm = TRUE)
	gcmrc3$bspec$CountF <- rowSums(gcmrc3$bspec[,40:56], na.rm = TRUE)
	gcmrc3$bspec$CountTotal <- rowSums(gcmrc3$bspec[,c(6:37, 40:56)], na.rm = TRUE)
	colnames(gcmrc3$bspec)[which(colnames(gcmrc3$bspec)=='Cpt5')] <- 'C0'
	colnames(gcmrc3$bspec)[which(colnames(gcmrc3$bspec)=='Fpt5')] <- 'F0'
	ccols <- paste0('C', 0:30)
	fcols <- paste0('F', 0:15)
	gcmrc4 <- gcmrc3
	gcmrc4$bspec <- gcmrc3$bspec[, c('BarcodeID', 'SpeciesID', ccols, 'CExtra', fcols, 'FExtra', 
		'CountC', 'CountF', 'CountTotal', 'Notes')]
	## Convert NA counts to 0s
	gcmrc4$bspec[is.na(gcmrc4$bspec)] <- 0
		gcmrc4$bspec <- droplevels(gcmrc4$bspec)	
	## Merge benthic sample and specimen data together. Same for drift
		d1 <- merge(gcmrc4$dsamp, gcmrc4$dspec, by = 'BarcodeID')
		b1 <- merge(gcmrc4$bsamp, gcmrc4$bspec, by = 'BarcodeID')	
	## Get densities and concentrations
	d1$Concentration <- d1$CountTotal / d1$Volume
	b1$Density <- b1$CountTotal / b1$SampleArea
	## Write data
	write.csv(d1, 'Data/DriftData.csv', row.names = FALSE)
	write.csv(b1, 'Data/BenthicData.csv', row.names = FALSE)
	write.csv(gcmrc4$sppl, 'Data/SpeciesList.csv', row.names = FALSE)

}
#BAreadwrite()

## Get data from GitHub
gitdat <- 'https://raw.githubusercontent.com/jmuehlbauer-usgs/BrightAngel/master/Data/'
w1 <- read.csv(paste0(gitdat, 'WhitingData.csv'))
wspp1 <- read.csv(paste0(gitdat, 'WhitingSpeciesList.csv'))
d1 <- read.csv(paste0(gitdat, 'DriftData.csv'))
b1 <- read.csv(paste0(gitdat, 'BenthicData.csv'))
spp1 <- read.csv(paste0(gitdat, 'SpeciesList.csv'))


##### Clean up data #####

## Convert dates and times to usable formats
d1$Date <- as.Date(d1$Date, format = '%Y-%m-%d')
	d1$ProcessDate <- as.Date(d1$ProcessDate, format = '%m/%d/%Y')
b1$Date <- as.Date(b1$Date, format = '%Y-%m-%d')
	b1$ProcessDate <- as.Date(b1$ProcessDate, format = '%Y-%m-%d')
w1$Date <- as.Date(w1$SampleDate, format = '%m/%d/%Y')	

## Add functional feeding groups to all datasets
d1$FFG <- spp1[match(d1$SpeciesID, spp1$SpeciesID), 'FFG']
b1$FFG <- spp1[match(b1$SpeciesID, spp1$SpeciesID), 'FFG']
w1$FFG <- spp1[match(w1$SpeciesID, spp1$SpeciesID), 'FFG']


##### Group specimens by FFG #####

## Combine data by FFG for benthics, drift, and Whiting
dffg <- tapply(d1$Concentration, d1$FFG, function(x){sum(x, na.rm = TRUE)})
bffg <- tapply(b1$Density, b1$FFG, function(x){sum(x, na.rm = TRUE)})
wffg <- tapply(w1$Density, w1$FFG, function(x){sum(x, na.rm = TRUE)})

## Combine all data streams, sort by trophic level
ffg <- as.data.frame(cbind(dffg, bffg, wffg))
	rownames(ffg)[1] <- 'Unknown'
ffg1 <- ffg[rownames(ffg) %in% c('Shredder', 'CollectorFilterer', 'CollectorGatherer', 'ScraperGrazer', 'Generalist', 'Predator'),]
ffg2 <- ffg1[c(6, 1, 2, 5, 3, 4),]
	colnames(ffg2) <- c('Drift', 'Benthic', 'Whiting')



### Stopped here. What follows is unverified and probably needs fixing. Above need to combine life stages for d1, b1, w1.

panel <- function(){}
par(mfrow = c(2, 1), mar = c(1.5, 5.5, 0.1, 0.1), oma = c(3.2, 0, 0, 0), xpd = FALSE)


wbar <- barplot(ffg2$Whiting, axes = FALSE, xlab = '', ylab = '', names.arg = FALSE, col = 'tomato')
box(bty = 'l')
axis(1, at = wbar, labels = rep('', length(wbar)))
axis(2, las = 2)
mtext(side = 2, expression(paste('Density (# * ', m^-2, ')')), line = 4)
legend('topleft', legend = '2010-2011', bty = 'n')


bbar <- barplot(ffg2$Benthic, axes = FALSE, xlab = '', ylab = '', names.arg = FALSE, col = 'slateblue1')
box(bty = 'l')
axis(1, at = bbar, labels = c('Shredders\n', 'Collector-\nfilterers', 'Collector-\ngatherers', 'Scrapers/\nGrazers', 'Generalists\n', 'Predators\n'), padj = 0.5)
axis(2, las = 2)
mtext(side = 2, expression(paste('Density (# * ', m^-2, ')')), line = 4)
mtext(side = 1, 'Functional feeding group', line = 3.5)
legend('topleft', legend = '2016-2017', bty = 'n')




###Calculating relative densities

ffg2$BenthicRel <- round(ffg2$Benthic / sum(ffg2$Benthic), 4)
ffg2$WhitingRel <- round(ffg2$Whiting / sum(ffg2$Whiting), 4)
#Note that I've combined a couple functions onto a single line of code (rather than calculating the sum on a separate line as I suggested over the phone). The result is the same, this is just a little cleaner. I've also included the round function, which just rounds the results the the specified number of decimal places (4 in this case).

#Subtracting Whitings from ours to see if there was an increase or decrease in densities.
ffg2$BenthicRel - ffg2$WhitingRel
Reldif <- (ffg2$BenthicRel - ffg2$WhitingRel) 


###Barplotting the difference in relative densities 
Relbar <- barplot(Reldif, axes = FALSE, xlab = '', ylab = '', names.arg = FALSE, col = 'slateblue1')
box(bty = 'l')
axis(1, at = bbar, labels = c('Shredders\n', 'Collector-\nfilterers', 'Collector-\ngatherers', 'Scrapers/\nGrazers', 'Generalists\n', 'Predators\n'), padj = 0.5)
axis(2, -1:1, las = 2)
mtext(side = 2, expression(paste('Density differences')), line = 4)
mtext(side = 1, 'Functional Feeding Groups', line = 3.5)
legend('topleft', legend = '', bty = 'n')


###ANOVA test
aov(abundance~ffg)
#error code. Begin here next time. 
 

