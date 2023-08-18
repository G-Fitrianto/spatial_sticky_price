setwd("D:/HARGA PASAR SPASIAL")
library(sp); library(sf); library(spdep); library(geosphere)

data <- read.csv("D:/HARGA PASAR SPASIAL/Harga Komoditas DIY Jateng Jatim Geocode.csv")

colnames(data)
# [1] "TANGGAL"              "PROVINSI"             "KOTA"                
# [4] "PASAR"                "KELOMPOK"             "JENIS"               
# [7] "HARGA"                "ALAMAT.LENGKAP.PASAR" "ALAMAT.PASAR"        
#[10] "LATITUDE"             "LONGITUDE"

data <- as.matrix(data)
class(data)
# [1] "matrix" "array"

data1        <- data.frame(rbind(matrix(c(rep(0, ncol(data))), ncol = ncol(data)),data))
pasar.id     <- which(data1$PASAR[2:nrow(data1)] != data1$PASAR[2:nrow(data1)-1], arr.ind = T)
pasar.m      <- matrix(NA, nrow = length(pasar.id), ncol = 2)
pasar.m[,1]  <- as.numeric(data1$LONGITUDE[pasar.id])
pasar.m[,2]  <- as.numeric(data1$LATITUDE[pasar.id])

# Distance and Inverse Distance Spatial W Matrix

dist.pasar     <- matrix(NA, length(pasar.id), length(pasar.id))
inv.dist.pasar <- matrix(NA, length(pasar.id), length(pasar.id))

for (j1 in 1:length(pasar.id)) {
  for (j2 in 1:length(pasar.id)) {
    dist.pasar[j1,j2]     <- distHaversine(pasar.m[j1,],pasar.m[j2,])/1000
    inv.dist.pasar[j1,j2] <- 1/(distHaversine(pasar.m[j1,],pasar.m[j2,]))*1000
  }
}

inf.detect <- which(inv.dist.pasar == Inf, arr.ind = T)
inv.dist.pasar[inf.detect] <- 0
# ----------------------------- #


# Nearest Neighbour Spatial W Matrix

k_neigh <- 5

points <- SpatialPoints(pasar.m)

nb.pasar.knn <- knn2nb(knearneigh(points, k_neigh))
W.pasar.knn  <- nb2mat(nb.pasar.knn, zero.policy = T, style = "W")
B.pasar.knn  <- nb2mat(nb.pasar.knn, zero.policy = T, style = "B")
# ----------------------------- #

date.det  <- which(data1$TANGGAL == data1$TANGGAL[2], arr.ind = T)
date.id   <- data1$TANGGAL[date.det[1]:(date.det[2]-1)]

# ----------------------------- #

comm.det      <- which(data1$JENIS[2:nrow(data1)] != data1$JENIS[2:nrow(data1)-1], arr.ind = T)
comm.id       <- data1$JENIS[comm.det]
comm.id       <- comm.id[-1]
comm.cat.ind  <- which(comm.id == comm.id[1], arr.ind = T)
comm.cat      <- comm.id[comm.cat.ind[1]:comm.cat.ind[2]]

# ----------------------------- #
# Berdasarkan Jenis Barang      #
# Agregasi Data by Time         #
# Analisis CGS                  #


comm.var  <- numeric(length(comm.cat))
path0     <- "D:/HARGA PASAR SPASIAL/Best Models by Commodity/"

delta        <- seq(0.001,0.999, 0.001)
W_mat        <- inv.dist.pasar

AIC_tab   <- matrix(NA, nrow = length(delta), ncol = length(comm.cat))

for (p1 in 1:length(comm.cat)) {
  data1.red         <- data1[-1,]
  data1.red.det     <- which(data1.red$JENIS == comm.cat[p1], arr.ind = T)
  data1.J1          <- data1.red[data1.red.det,]
  data1.date.det    <- which(data1.J1$TANGGAL == date.id[1], arr.ind = T)
  data1.J1.date.red <- data1.J1[data1.date.det,]
  
  ## Estimasi ##
  
  X_data      <- matrix(NA, ncol = 1, nrow = length(pasar.id))
  date_length <- length(date.id) 
  
  for (m1 in 1:length(pasar.id)) {
    min_id = ((m1-1)*date_length) + 1
    max_id = m1*date_length
    
    X_data[m1] <- mean(na.omit(as.numeric(data1.J1$HARGA[min_id:max_id])))
  }
  
  nan.det  <- which(X_data == "NaN", arr.ind = T)
  Y_dot    <- matrix(X_data[-nan.det[,1]], ncol = 1)
  
  W_dot        <- W_mat[-nan.det[,1],-nan.det[,1]]
  I_dot        <- diag(1, nrow(Y_dot), nrow(Y_dot))
  i_dot        <- matrix(1, ncol = 1, nrow = nrow(Y_dot))
  WY_dot       <- W_dot %*% Y_dot
  delta_par    <- numeric(length = length(delta))

  countif    <- 0 
  for (k1 in 1:length(delta)) {
    d_val     <- delta[k1]
    B_mat     <- (I_dot - d_val*W_dot)
    Y_til     <- B_mat %*% Y_dot
    WY_til    <- B_mat %*% i_dot
    
    model1          <- lm(Y_til~WY_til)
    delta_par[k1]   <- coef(model1)[2]
    AIC_tab[k1,p1]  <- AIC(model1)
  }
  
  # min_eigenW <- abs(min(eigen(W_dot)$values))
  # max_eigenW <- abs(max(eigen(W_dot)$values))
  # plot(AIC_par, pch=16, cex=2, col="red")#, ylim=c(min(AIC_par), max(AIC_par)))
  
  minAIC       <- which(AIC_tab[,p1] == min(AIC_tab[,p1]), arr.ind = T)
  comm.var[p1] <- delta[minAIC]
  
  path2          <- paste0(path0, paste0(comm.cat[p1], paste0(" ","Estimation.txt")))
  d_val.best     <- delta[minAIC]
  B_mat.best     <- (I_dot - d_val.best*W_dot)
  Y_til.best     <- B_mat.best %*% Y_dot
  WY_til.best    <- B_mat.best %*% i_dot
  
  model1.best    <- lm(Y_til.best~WY_til.best)
  best.summary   <- capture.output(summary(model1.best))
  
  delta.best     <- paste0("Spatial Error Parameter = ", comm.var[p1])
  writeLines(c(delta.best, best.summary), path2)
}

# ------------------------------------- #
# Berdasarkan Jenis Barang dan Time     #


start.time <- Sys.time()

comm.var2 <- numeric(length(comm.cat))
path0.1   <- "D:/HARGA PASAR SPASIAL/Best Models by Commodity and Time/"

delta     <- seq(0.01,0.99, 0.01)
W_mat     <- inv.dist.pasar
I_mat     <- diag(1, length(date.id), length(date.id))
i_vec     <- matrix(1, nrow = length(date.id), ncol = 1)
A_mat     <- kronecker(I_mat, W_mat)

AIC_tab2  <- matrix(NA, nrow = length(delta), ncol = length(comm.cat))


for (p2 in 1:length(comm.cat)) {
  data1.red         <- data1[-1,]
  data1.red.det     <- which(data1.red$JENIS == comm.cat[p2], arr.ind = T)
  data1.J1          <- data1.red[data1.red.det,]
  
  ## Estimasi ##
  
  XT_data      <- matrix(as.numeric(data1.J1$HARGA), ncol = 1)
  
  na.det    <- which(is.na(XT_data), arr.ind = T)
  XT_dot    <- matrix(XT_data[-na.det[,1]], ncol = 1)
  
  A_dot        <- A_mat[-na.det[,1],-na.det[,1]]
  I_dot        <- diag(1, nrow(XT_dot), nrow(XT_dot))
  i_dot        <- matrix(1, ncol = 1, nrow = nrow(XT_dot))
  WY_dot       <- A_dot %*% XT_dot
  delta_par2   <- numeric(length = length(delta))
  
  countif    <- 0 
  for (k1 in 1:length(delta)) {
    d_val     <- delta[k1]
    B_mat     <- (I_dot - d_val*A_dot)
    XT_til    <- B_mat %*% XT_dot
    WY_til    <- B_mat %*% i_dot
    
    model2           <- lm(XT_til~WY_til)
    delta_par2[k1]    <- coef(model2)[2]
    AIC_tab2[k1,p2]  <- AIC(model2)
  }
  
  # min_eigenW <- abs(min(eigen(W_dot)$values))
  # max_eigenW <- abs(max(eigen(W_dot)$values))
  # plot(AIC_par, pch=16, cex=2, col="red")#, ylim=c(min(AIC_par), max(AIC_par)))
  
  minAIC       <- which(AIC_tab2[,p2] == min(AIC_tab2[,p2]), arr.ind = T)
  comm.var2[p2] <- delta[minAIC]
  
  path2.1        <- paste0(path0.1, paste0(comm.cat[p2], paste0(" ","Estimation.txt")))
  d_val.best     <- delta[minAIC]
  B_mat.best     <- (I_dot - d_val.best*A_dot)
  XT_til.best    <- B_mat.best %*% XT_dot
  WY_til.best    <- B_mat.best %*% i_dot
  
  model2.best    <- lm(XT_til.best~WY_til.best)
  best.summary   <- capture.output(summary(model2.best))
  
  delta.best     <- paste0("Spatial Error Parameter = ", comm.var2[p2])
  writeLines(c(delta.best, best.summary), path2.1)
}

end.time <- Sys.time()
elapsed  <- end.time - start.time
# ------------------------------------- #