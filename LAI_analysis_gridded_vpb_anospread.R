##############################################################################
# converts non-gridded L2 GHG data with grid infromation to gridded GHG data #
# a temporal averaging of the data of one grid cell is performed             #
##############################################################################

#############
# LIBRARIES #
#############

#library(cmsaf)
library(ncdf4)
library(RNetCDF)
library(fields)  
library(rgdal)
library(sp)

#library(maptools)
#library(RColorBrewer)
#library(colorRamps)
#require(mapdata)
#library(rworldmap) 
#library(rworldxtra)  
#library(ggplot2)

###############
# Definitions #
###############

# INPUT
indir <- "/cmsaf/cmsaf-rad6/lkeupp/politik/LSA-SAF"
#inname <- "LSASAF_MSG_LAI_Euro_2014-2015_grid_time.nc"
#inname <- "LSASAF_MSG_LAI_Euro_2014-2015_gridtest.nc"
#inname <- "LSASAF_MSG_LAI_Euro_2014-2015_interp_1deg_masked.nc"
inname <- "LSASAF_MSG_LAI_Euro_2006-2015_interp_masked.nc"
years <- c(2006:2015)
cty <- length(years)
infile <- file.path(indir,inname)
#infile2 <- file.path(indir,inname2)

var <- "LAI" 
un <- "m^2*(m^-2)"
latname <- "lat"
lonname <- "lon"

# for regional mean time series
region <- TRUE
lomin <- 6
lomax <- 15
lamin <- 47
lamax <- 55
onamer <- "LAI_regional"


# OUTPUT
#outdir <- "/cmsaf/cmsaf-rad6/lkeupp/politik/LSA-SAF/out_2006-2015_n1d"
outdir <- "/cmsaf/cmsaf-rad6/lkeupp/politik/LSA-SAF"

nd <- 1
ud <- (nd-1)/2


outname <- paste0("LSASAF_LAI_vpb_2006-2015_", nd,  "rm.nc")
outfile <- file.path(outdir,outname)





#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#
#---------------------------------------------------------------------------------------------------------------#
#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#


#########################################################################################
# read data, truncate (only one valid digit), make data.frame, save data.frame as table #
#########################################################################################
ncd <- nc_open(infile, write=FALSE, readunlim=FALSE, verbose=FALSE, auto_GMT=FALSE,suppress_dimvals=TRUE)
co <- ncvar_get(ncd, varid=var, start=NA, count=NA, verbose=FALSE, signedbyte=TRUE, collapse_degen=TRUE)
tim <- ncvar_get(ncd, varid='time', start=NA, count=NA, verbose=FALSE, signedbyte=TRUE, collapse_degen=TRUE)
#un <- ncatt_get(ncd, var, attname='units', verbose=FALSE)
scf <- ncatt_get(ncd, var, attname='scale_factor', verbose=FALSE)
la <- ncvar_get(ncd, varid=latname, start=NA, count=NA, verbose=FALSE, signedbyte=TRUE, collapse_degen=TRUE)
lo <- ncvar_get(ncd, varid=lonname, start=NA, count=NA, verbose=FALSE, signedbyte=TRUE, collapse_degen=TRUE)

nlo <- length(lo)
nla <- length(la)
#londim <- ncd$dim[["lon"]]
#lon <- londim$vals
#nlo <- londim$len
#latdim <- ncd$dim[["lat"]]
#lat <- latdim$vals
#nla <- latdim$len

timedim <- ncd$dim[["time"]] 
nt <- timedim$len
time.unit <- timedim$units
time <- timedim$vals

date.time <- as.Date(utcal.nc(time.unit,tim,type="s"))
#date.time2 <- as.Date(utcal.nc(time.unit,tim,type="n"))
#timeday <- format(date.time,"%d")
#timemonth <- format(date.time,"%m")
#timeyear <- format(date.time,"%Y") 

nc_close(ncd)

print(paste(scf,sep=" "))
print(paste(co[1]))
print(paste(la[1]))
print(paste(lo[1]))

print(paste("Data read"))

##########################################


#sc <- scf$value
#co <- co/(sc*sc)
co[co < 0] <- NA

mv <- -9e+06

##########################################################

#######################
# alternative dates
#######################

# 2015.05451
#year.month(day as portion of month)

teti <- utcal.nc(time.unit,tim)
teti <- teti[,1:3]
datecom <-  teti[,3]
datecom[] <- mv
ty <- c(1:12)
ty <- rep(ty, cty)
ty <- matrix(ty, ncol=cty)
daynum <- ty
daynum[,] <- 0

	for(y in 1:cty) {
	ye <- 2006+y-1
	mon <- 0


	for( m in 1:12) {
	mon <- mon+1


for(t in 1:nt) {

    if(teti[t,1] == ye){
    if(teti[t,2] == mon){

year <- as.character(teti[t,1])
daynum[m,y] <- daynum[m,y]+1


mo <- teti[t,2]
if(mo < 10) {
month <- paste0("0", as.character(mo))
} else {
month <-  as.character(mo)
}

ty[m,y] <- paste(year, month, sep=".")


    }
    }


}
	}
	}

for(t in 1:nt) {

y <- (teti[t,1]-2006+1)
m <- teti[t,2]
d <- teti[t,3]

if(y == 1 && m == 1) {
day <- as.numeric((d-1)/31)
}else {
day <- as.numeric((d-1)/daynum[m,y])
}

trunc1 <- function(x, digits) trunc(x*10^digits)/10^digits
#print(paste(t, y, m, d, day, trunc1(day*1000,0), sep=", "))
day <- trunc1(day*1000,0)

if(day < 100) {

day <- paste0("0", as.character(day))

} else {

day <- as.character(day)

}
datecom[t] <- paste0(ty[m,y], day)
datecom[t] <- as.numeric(datecom[t]) 

}



# year.(days since beginning of the year)
teti <- utcal.nc(time.unit,tim)
teti <- teti[,1:3]  

daynu <-  c(1:cty)
dayn <- c(1:nt)
daynu[] <- 0
dayn[] <- 0

        for(y in 1:cty) {
        ye <- 2006+y-1
	dcy <- 0

for(t in 1:nt) {

    if(teti[t,1] == ye){ 
dcy <- dcy+1

year <- as.character(teti[t,1])
daynu[y] <- dcy
dayn[t] <- dcy

    }
   
 
 
}
        }
     

datecom2 <- (teti[,1]*1000)+dayn

for(t in 1:nt) {

y <- (teti[t,1]-2006+1)
if(y == 1) {
datecom2[t] <- datecom2[t]+3
}


}

datecom2 <- datecom2/1000

#############################################################################################

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Data inside a region
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
if(region == TRUE) {
#-#-#-#-#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#--#-#
loct <- 0
lact <- 0

for(o in 1:nlo){
	if(lo[o] < lomin){
#		print(paste(lo[o], "<", lomin, sep=" "))
	} else if (lo[o] >= lomin && lo[o] <= lomax) {
		loct <- loct+1
		if(loct == 1){
			xb <- o
#			print(paste(lo[xb], "=", lomin, sep=" "))
		} else {
			xe <- o
#			print(paste(lo[xe], "=", lomax, sep=" "))
		}
	} else {
#		print(paste(lo[o], ">", lomax, sep=" "))
	}
}

print(paste("lo[", xb, "] = ",lo[xb], ", lomin = ", lomin, sep=" "))
print(paste("lo[", xe, "] = ",lo[xe], ", lomax = ", lomax, sep=" "))



for(o in 1:nla){
	if(la[o] < lamin){
#		print(paste(la[o], "<", lamin, sep=" "))
	} else if (la[o] >= lamin && la[o] <= lamax) {
		lact <- lact+1
		if(lact == 1){
			yb <- o
#			print(paste(la[yb], "=", lamin, sep=" "))
		} else {
			ye <- o
#			print(paste(la[ye], "=", lamax, sep=" "))
		}
	} else {
#		print(paste(la[o], ">", lamax, sep=" "))
	}
}

print(paste("la[", yb, "] = ",la[yb], ", lamin = ", lamin, sep=" "))
print(paste("la[", ye, "] = ",la[ye], ", lamax = ", lamax, sep=" "))


cor <- co[c(xb:xe),c(yb:ye),]

print(paste("Regional delineated"))
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
} else {
cor <- co
xb <- 1
xe <- length(cor[,1,1])
yb <- 1
ye <- length(cor[1,,1])

}
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#####################################################################

xa <- 1
xz <- length(cor[,1,1])
ya <- 1
yz <- length(cor[1,,1])

ngc <- xz*yz*cty

vpb <- rep(mv, ngc)
vpb <- array(vpb, dim=c(xz,yz,cty))
vpbn <- vpb

for(xct in xa:xz) {		#x
for(yct in ya:yz) {		#y

#for(xct in 4:4) {
#for(yct in 8:8) {



print(paste(xct, lo[xct], yct, la[yct], sep=", "))


#########################################################
# running mean per grid cell
##########################################################

cou <- cor[xct,yct,]
	ctna <- sum(is.na(cou))
	if(ctna == nt) {
		vpb[xct,yct,] <-mv
		print(paste0("no data for: lon ", lo[xct], " - lat ", la[yct]))
	} else {

courm <- cou
#cou3rm <- mv

beg <- 1+ud
fin <- nt-ud
courm[1:ud] <- mv

if(fin != nt) {
fint <- fin+1
courm[fint:nt] <- mv
}


for (i in beg:fin) {
	st <- i-ud
	en <- i+ud
	ndna <- nd-sum(is.na(cou[c(st:en)]))
	courm[i] <- sum(cou[c(st:en)], na.rm=T)/ndna

if(is.nan(courm[i]) == TRUE) { is.na(courm[i]) <- TRUE }

}

courm[courm==mv] <- NA
maxabs <- max(courm, na.rm=TRUE)
minabs <- min(courm, na.rm=TRUE)

###########################################################
# find absolute minimum per year
###########################################################

absmins <- rep(0, cty)

for(z in 1:cty) {

a <- 1
b <- 150 #daynu[z]
if(z == 1){
b <- b-3
}

z1 <- z-1
if(z > 1){
a <- a+sum(daynu[1:z1])
b <- b+sum(daynu[1:z1])
#a <- a+(z-1)*365
#b <- b+(z-1)*365
}

amin <- courm[a:b]
amin[] <- 0
amv <- min(courm[a:b], na.rm=T)

for (m in a:b) {
	if (is.na(courm[m]) == FALSE && is.nan(courm[m]) == FALSE) {
	if(courm[m] == amv){
		print(paste("absolute minimum at", m, "(", date.time[m], ") :", courm[m], amv))
		amin[m] <- m
		}
	}
}

	

amin[is.na(amin)] <- 0

if( length(amin[amin!=0]) > 1) {
amino <- amin[amin!=0]
absmins[z] <- amino[length(amino)]
print(paste('more than on absolute minimum', z, xct, yct, amino, courm[amino], sep=", "))
} else if(sum(amin)==0) {
absmins[z] <- NA 
} else {
absmins[z] <- amin[amin!=0]
#amin <- amin[amin!=0]
}
}



vpbn[xct,yct,] <- absmins
vpb[xct,yct,] <- (datecom2[absmins]-years)*1000

}  	#if(countna !0 = nt)
}	#x
}	#y



################################################################################################


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#timy <- date.time
timy <- datecom2
#timy <- datecom
#timy <- 2000+tim/365
#timy <- tim
############################################################################


#vpb <- (datecom2[vpb]-years)*1000

#time <- mean(tim)
#time <- tim #[n]
time <- years-2000
time.unit <- "years since 2000-01-01"
un <- "day of year"
var <- "begin of vegetation period from LAI"
lon <- lo[c(xb:xe)] 
lat <- la[c(yb:ye)]
vpb <- vpb


x <- ncdim_def(name="lon",units="degrees_east", vals=lon)
y <- ncdim_def(name="lat",units="degrees_north", vals=lat)
t <- ncdim_def(name="time",units=time.unit, vals=time, unlim=FALSE) #TRUE)
var1 <- ncvar_def(var,as.character(un),list(x,y,t),mv,prec="float")
#var1 <- ncvar_def(var,as.character(un$value),list(x,y,t),mv,prec="float")
vars <- list(var1)   
ncnew <- nc_create(outfile,vars)
ncvar_put(ncnew,var1,vpb)
ncatt_put(ncnew,"lon","standard_name","longitude",prec="text")
ncatt_put(ncnew,"lat","standard_name","latitude",prec="text")
nc_close(ncnew)

#ncinfo(outfile,"m")

print(paste('Netcdf file written'))



vpbmx <- vpb[,,1]
vpbmx[,] <- mv
vpbmn <- vpbmx
 
for(xct in xa:xz) { 
for(yct in ya:yz) {

vpbmx[xct,yct] <- max(vpb[xct,yct,])
vpbmn[xct,yct] <- min(vpb[xct,yct,])

}
}


#nacor <- rep(mv, nt)
#tpn <- 0
#
#for(i in 1:nt){
#nacor[i] <- sum(is.na(cor[,,i]))
#if(nacor[i] > 12){
#print(i)
#tpn <- c(tpn, i)
#}
#
#}
#tpn <- tpn[2:length(tpn)]


min(cor[4,8,sum(daynu[1:7]):sum(daynu[1:8])], na.rm=T)

for(yct in 1:8) {
png(file=paste0(outdir,"/vpb_lai_", as.character(nd), "rm", as.character(yct), ".png"), width=960)

colo <- c("red", "green", "orange", "blue", "purple", "darkgreen", "magenta", "cyan", "gray")
tit <- paste0("Zeitreihe des Beginns der Vegetationsperiode berechnet aus LAI")
sub <- paste0("y= ", yct, " (", lat[yct], "°N)")
plot(2006:2015, vpb[1,yct,], type="b", ylim=c(0,150), col=colo[1], main=tit, sub=sub)
lines(2006:2015, vpb[2,yct,], type="b", col=colo[2])
lines(2006:2015, vpb[3,yct,], type="b", col=colo[3])
lines(2006:2015, vpb[4,yct,], type="b", col=colo[4])
lines(2006:2015, vpb[5,yct,], type="b", col=colo[5])
lines(2006:2015, vpb[6,yct,], type="b", col=colo[6])
lines(2006:2015, vpb[7,yct,], type="b", col=colo[7])
lines(2006:2015, vpb[8,yct,], type="b", col=colo[8])
lines(2006:2015, vpb[9,yct,], type="b", col=colo[9])
text(x=2006, y=0, lab=trunc1(lon[1],2), col=colo[1], pos=4)
text(x=2007, y=0, lab=trunc1(lon[2],2), col=colo[2], pos=4)
text(x=2008, y=0, lab=trunc1(lon[3],2), col=colo[3], pos=4)
text(x=2009, y=0, lab=trunc1(lon[4],2), col=colo[4], pos=4)
text(x=2010, y=0, lab=trunc1(lon[5],2), col=colo[5], pos=4)
text(x=2011, y=0, lab=trunc1(lon[6],2), col=colo[6], pos=4)
text(x=2012, y=0, lab=trunc1(lon[7],2), col=colo[7], pos=4)
text(x=2013, y=0, lab=trunc1(lon[8],2), col=colo[8], pos=4)
text(x=2014, y=0, lab=trunc1(lon[9],2), col=colo[9], pos=4)
text(x=2014.75, y=0, lab="°E", col="black", pos=4)
dev.off()

}



#######
# END #
#######

d
