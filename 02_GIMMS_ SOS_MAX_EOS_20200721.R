# Calculation SOS, MAX, EOS adjusting growing season by centering MAX
# Author: Roberto O. Chávez
# Last update: 23-06-2020

library(raster)
library(npphen)
library(rts)
library(lubridate)

# inpath, outpath
inpath <- "~/PROJECTS/01b_phen_chile_GIMMSv1/01_v1_ndvi_tif"
ndvi.list <- list.files(path=inpath, pattern=glob2rx("Chile_gimms*.tif"), full.names=T)
dates.table <- read.csv("~/PROJECTS/01b_phen_chile_GIMMSv1/table/GIMMS_NDVI_3gv1_828_dates_35_GS.csv", sep = ",", header=TRUE)
# dates
dd<-as.character(dates.table$dd)
mm<-as.character(dates.table$mm)
YY<-as.character(dates.table$YY)
ddmmYY<-paste(dd,"/",mm,"/",YY, sep="")
all_dates <- as.Date(ddmmYY, format='%d/%m/%Y')

#======================================================================================================================================
# Extracting a sample ts
vi.st <- stack(ndvi.list)
vi.rts<-rts(vi.st, all_dates)
px<-cellFromXY(vi.rts,c( -71.36,-37.85)) 
px2<-cellFromXY(vi.rts,c( -71.370,-31.121))

vi.ts<-extract(vi.rts,px)
data<-as.data.frame(vi.ts)
data$dmy <- as.Date(all_dates, format='%d-%m-%Y')
ndvi <- as.numeric(data$V1)
data$ndvi <- ndvi

vi.ts2<-extract(vi.rts,px2)
data2<-as.data.frame(vi.ts2)
data2$dmy <- as.Date(all_dates, format='%d-%m-%Y')
ndvi <- as.numeric(data2$V1)
data2$ndvi <- ndvi

setwd("~/PROJECTS/01b_phen_chile_GIMMSv1/table")
write.csv(data,"px1.csv")
write.csv(data2,"px2.csv")

data <- read.csv("px1.csv")
data$dmy <- as.Date(all_dates, format='%d-%m-%Y')

data2 <- read.csv("px2.csv")
data2$dmy <- as.Date(all_dates, format='%d-%m-%Y')
#----------------------------------------------------------------------------------------

# Function for a single vector

x=data2$ndvi
dates=data2$dmy
rge=c(0,1)

Phen_gs_adj_map <- function (x, dates, rge, th) { ### NO DEBERIAMOS USAR nGS EN ESTA FUNCION, YA LO SAQUE....
  #h=2 # setea el hemisferio 2 # No conviene tener variables con valores fijos...
  options(warn=-1) # elimina los warnings de la creación de la tabla de DOY y DGS
  if (length(rge) != 2) {
    stop("rge must be a vector of length 2")
  }
  if (rge[1] > rge[2]) {
    stop("rge vector order must be minimum/maximum")
  }
  if (length(dates) != length(x)) {
    stop("N of dates and files do not match")
  }
  if (all(is.na(x))) {          ### MATI, ESTA PARTE ES CLAVE, CUAL ES LA LONGITUD DEL OUTPUT VECTOR????
    return(rep(NA, 8))
  }
  DOY <- yday(dates)
  D1 <- cbind(DOY, x)
  if (length(unique(D1[, 2])) < 10 | (nrow(D1) - sum(is.na(D1))) < (0.1 * length(D1))) {
    return(rep(NA, 8))        ### MATI AQUI LO MISMO, ESTA PARTE ES CLAVE, CUAL ES LA LONGITUD DEL OUTPUT VECTOR????
  }
  DOGS <- cbind(seq(1, 365), c(seq(185, 365), seq(1, 184)))     ##### AQUI LO OBLIGAMOS A SER H=2 (HEMISFERIO SUR) SIEMPRE
  for (i in 1:nrow(D1)) {
    D1[i, 1] <- DOGS[which(DOGS[, 1] == D1[i, 1], arr.ind = TRUE),2]
  }
  Hmat <- Hpi(na.omit(D1))
  Hmat[1, 2] <- Hmat[2, 1]
  K1 <- kde(na.omit(D1), H = Hmat, xmin = c(1, rge[1]), xmax = c(365,rge[2]), gridsize = c(365, 500))
  K1Con <- K1$estimate
  for (j in 1:365) {
    K1Con[j, ] <- K1$estimate[j, ]/sum(K1$estimate[j, ])
  }
  MAXY <- apply(K1Con, 1, max)
  for (i in 1:365) {
    MAXY[i] <- median(K1$eval.points[[2]][which(K1Con[i,] == MAXY[i], arr.ind = TRUE)])
  }
  #first.DOGS <- which.min(D1[, 1])
  #first.DOGS <- first.DOGS[1]
  #last.DOGS <- first.DOGS + nGS - 1
  Ref <- rep(NA, 365)
  for (i in 1:365) {
    Ref[i] <- MAXY[i]
  }
  dGS <- seq(1, 365, 1) ###### RCH NUEVO
  
  # Plotting
  plot(dGS, Ref[1:365], 
       xlab = "DGS", ylab = "VI", font.lab = 2, type = "l")
  
  #---------------------- Checking max -------------------------------------
  
  fila <- Ref[1:365] ###### MO NUEVO = CREA EN UN VECTOR LOS DATOS REORDENADOS
  
  # Checking if there are more than one max
  p <- max(fila) # Calcula el/los máximo(s)
  
  ## NO ENTENDI ESTA PARTE, ASI QUE LA DESACTIVE
  #pp <- floor(length(fila[which(fila==p)])/2)+1 # Crea un vector de extensión 50%+1 para centrar
  
  
  #if (pp>1){ # Si hay más de un máximo, deje el DOY del centro (50%+1)
  #  peak.vector <- which.max(fila)+pp-1
  #  dGS_max <- dGS[peak.vector]
  #} else { # Si no, deje el DOY del único valor máximo
  #  peak.vector <- which.max(fila)
  #  dGS_max <- dGS[peak.vector]
}
# Plotting
# abline(v=dGS_max, col="red") # Marca el DOY del peak

dGS_max <- dGS[which.max(Ref[1:365])]
abline(v=dGS_max, col="red")
  
  
  # Reordena las fechas para central al peak
  le_A <- 365-dGS_max+1 
  ini_A <- 182
  end_A <- 182+le_A-1
  if(end_A>365){
    sobra <- end_A
    end_A <- 365
    ini_B <- 1
    le_B <- sobra-365
    le_A <- le_A-le_B
    end_B <- ini_B+le_B-1
    if((le_A+le_B)<365) {
      ini_C <- end_B+1
      end_C <- 181
      SQ <- c(seq(ini_C,end_C),seq(ini_A,end_A),seq(ini_B,end_B))
    }
  } else {
    ini_B <- end_A + 1
    end_B <- 365
    le_B <- end_B - ini_B + 1
    if((le_A+le_B)<365) {
      le_C <- 365-le_A-le_B-1
      ini_C <- 1
      end_C <- le_C
      SQ <- c(seq(ini_B,end_B),seq(ini_C,end_C),seq(ini_A,end_A))
    }
  }
  
  ######### HASTA AQUI NO TOCAR
  newDOGS <- cbind(seq(1, 365), SQ) #Esta es la linea que arroja warnings
  
  for (i in 1:nrow(D1)) {
    D1[i, 1] <- newDOGS[which(newDOGS[, 1] == D1[i, 1], arr.ind = TRUE),2]
  }
  Hmat <- Hpi(na.omit(D1))
  Hmat[1, 2] <- Hmat[2, 1]
  K1 <- kde(na.omit(D1), H = Hmat, xmin = c(1, rge[1]), xmax = c(365,rge[2]), gridsize = c(365, 500))
  K1Con <- K1$estimate
  for (j in 1:365) {
    K1Con[j, ] <- K1$estimate[j, ]/sum(K1$estimate[j, ])
  }
  MAXY <- apply(K1Con, 1, max)
  for (i in 1:365) {
    MAXY[i] <- median(K1$eval.points[[2]][which(K1Con[i,] == MAXY[i], arr.ind = TRUE)])
  }
  #first.newDOGS <- which.min(D1[, 1])
  #first.newDOGS <- first.newDOGS[1]
  #last.newDOGS <- first.newDOGS + nGS - 1
  Ref2 <- rep(NA, 365)
  for (i in 1:365) {
    Ref2[i] <- MAXY[i]
  }
  plot(seq(1, 365, 1), Ref2[1:365],xlab = "DGS", ylab = "VI", font.lab = 2, type = "l")
  Ref2[1:365]
  ndGS <- seq(1, 365, 1)
  ndGS_max <- ndGS[which.max(Ref2[1:365])] #Sacar
  abline(v=ndGS_max, col="blue")
  
  ####### HASTA AQUI ARREGLÉ #######################################################################
  
  # DOY-DGS1-DGS2 Conversion table
  c.table <- as.data.frame(cbind(DOGS[,2],DOGS[,1],newDOGS[,2]))#doy-dgs1-dgs2
  names(c.table)<- c("doy","dgs1","dgs2")
  
  
  # Calculating VI and DOY Phenometrics
  phen.adj <- Ref2[first.newDOGS:last.newDOGS]
  
  ## Checking if there are more than one max ... again
  set<- phen.adj[7:19] # Extrae un vector del centro de la phen, por si hay peaks separados, los hay.
  p2 <- max(set) # Calcula el peak del centro
  pp2 <- floor(length(phen.adj[which(set==p2)])/2) # Crea un vector de extensión 50%+1 para centrar
  
  if (pp2>1){ # Si hay más de un máximo, deje el DOY del centro (50%+1)
    peak.vector.2 <- which.max(set)+pp2-1
    peak.vector.3 <- 6 + peak.vector.2
    ndGS_max <- ndGS[peak.vector.3]
  } else { # Si no, deje el DOY del único valor máximo
    peak.vector.2 <- which.max(set)
    peak.vector.3 <- 6 + peak.vector.2
    ndGS_max <- ndGS[peak.vector.3]
  }
  
  # Phenometric: Peak
  doy.max <- c.table$doy[c.table$dgs2==ndGS_max]
  vi.max <- round(max(phen.adj)*10000)
  
  # Phenometric: SOS
  bef <- phen.adj[1:peak.vector.3] # Extrae un vector desde el inicio hasta el peak (incluido)
  sos.p <- ((max(phen.adj,na.rm = T)-min(bef,na.rm = T))*th)+min(bef,na.rm = T) # Calcula TH
  rbef <- rev(bef) # Invierte el orden del vector para buscar desde el peak 
  sos.val <- rbef[which(rbef<sos.p)-1][1] # Extrae el primer valor donde se cumpla < que el TH
  sos.tim <- last(which(rbef==sos.val)) # Extrae el primer valor en caso de que el TH sea mas de 1
  sos.time <- as.numeric(length(bef)-sos.tim+1) # Calcula el DGS en la serie original
  ndGS_sos <- ndGS[sos.time] # Extrae el DOY de la tabla
  
  doy.sos <-  c.table$doy[c.table$dgs2==ndGS_sos]
  vi.sos <- round(phen.adj[sos.time]*10000)
  
  if(vi.sos==vi.max){ # Si el VI.sos es igual al VI.Max
    ##
    bef <- phen.adj[1:peak.vector.3] # Extrae un vector desde el inicio hasta el peak (incluido)
    sos.p <- ((max(phen.adj,na.rm = T)-min(bef,na.rm = T))*th)+min(bef,na.rm = T) # Calcula TH
    rbef <- rev(bef) # Invierte el orden del vector para buscar desde el peak
    sos.val <- rbef[which(rbef<sos.p)-1][1] # Extrae el primer valor donde se cumpla < que el TH
    sos.tim <- last(which(rbef==sos.val)) # Extrae el primer valor en caso de que el TH sea mas de 1
    sos.time <- as.numeric(length(bef)-sos.tim) # Calcula el DGS anterior al peak
    ndGS_sos <- ndGS[sos.time] # Extrae el DOY de la tabla
    
    doy.sos <-  c.table$doy[c.table$dgs2==ndGS_sos]
    vi.sos <- round(phen.adj[sos.time]*10000)
    ###
  }
  
  # Phenometric: EOS
  aft <- phen.adj[peak.vector.3:nGS] # Extrae un vector desde el peak (incluido) hasta el final
  eos.p <- ((max(phen.adj,na.rm = T)-min(aft,na.rm = T))*th)+min(aft,na.rm = T) # Calcula TH
  eos.val <- aft[which(aft<=eos.p)[1]] # Extrae el primer valor donde se cumpla < que el TH
  eos.tim <- which(aft==eos.val)[1] # Busca el DGS en la serie corta
  eos.time <- (nGS-length(aft))+(as.numeric(eos.tim)) # Calcula el DGS en la serie original
  ndGS_eos <- ndGS[eos.time] # Extrae el DOY de la tabla
  
  doy.eos <-  c.table$doy[c.table$dgs2==ndGS_eos]
  vi.eos <- round(phen.adj[eos.time]*10000)
  
  if(vi.eos==vi.max){ # Si el VI.eos es igual al VI.Max
    aft <- phen.adj[peak.vector.3:nGS] # Extrae un vector desde el peak (incluido) hasta el final
    eos.p <- ((max(phen.adj,na.rm = T)-min(aft,na.rm = T))*th)+min(aft,na.rm = T) # Calcula TH
    eos.val <- aft[which(aft<=eos.p)[1]] # Extrae el primer valor donde se cumpla < que el TH
    eos.tim <- last(which(aft==eos.val)) # Busca el último DGS en la serie corta
    eos.time <- (nGS-length(aft))+(as.numeric(eos.tim)) # Calcula el DGS en la serie original
    ndGS_eos <- ndGS[eos.time] # Extrae el DOY de la tabla
    
    doy.eos <-  c.table$doy[c.table$dgs2==ndGS_eos]
    vi.eos <- round(phen.adj[eos.time]*10000)
    
  }
  
  #Plot check
  plot(y=phen.adj,x=dGS,type="b",main=z)
  abline(v=c(ndGS_sos,ndGS_max,ndGS_eos),col=c("blue","green","red"))
  abline(h=sos.p,lty=3,col="blue")
  abline(h=sos.val,col="blue")
  abline(h=eos.p,lty=3,col="red")
  abline(h=eos.val,col="red")
  
  # Phenometric: LGS
  lgs <- ndGS_eos-ndGS_sos
  
  # Phenometric: VIsum
  VIsum <- round(sum(phen.adj)*10000)
  
  ### Output 
  df.phen1px <- cbind.data.frame(doy.sos,doy.max,doy.eos,vi.sos,vi.max,vi.eos,lgs,VIsum)
  return(df.phen1px)
}

