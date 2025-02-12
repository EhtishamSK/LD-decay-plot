setwd("C:/Users/ehtis/OneDrive - New Mexico State University/SUNNY/Research Projects/Fruit Morphology Projects/GWAS/LD")
mydata <- read.table("LD_tassel_output.txt", header = TRUE)
str(mydata)
mydata$dist <- as.numeric(mydata$dist)
#creating a new column named rsq in the data frame data and populating it with the values from the existing column R2
myresult <- mydata
myresult$rsq <- mydata$R2
# Remove rows with NaN in the 'R2' column
file <- na.omit(myresult)
head(file)
# this data was obtained from genotyping of 128 individuals, so n = no. of gametes = no. of sampled chromosomes = 12 x 128
n = 1536
# C values range from about 0.5 to 2
Cstart <- c(C=0.5)

# let's fit a non linear model using the arbitrary C value
modelC <- nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))), data=file, start=Cstart, control=nls.control(maxiter=100))

# extract rho, the recombination parameter, 4Nr
rho <- summary(modelC)$parameters[1]

# feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
newrsq <- ((10+rho*file$dist)/((2+rho*file$dist)*(11+rho*file$dist)))*(1+((3+rho*file$dist)*(12+12*rho*file$dist+(rho*file$dist)^2))/(n*(2+rho*file$dist)*(11+rho*file$dist)))

newfile <- data.frame(file$dist, newrsq)

#maxld <- max(file$rsq) #using max LD value from initial input file
maxld <- max(newfile$newrsq) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecaydist <- newfile$file.dist[which.min(abs(newfile$newrsq-halfdecay))]
newfile <- newfile[order(newfile$file.dist),]

tiff("chr1_LD_DECAY4.tif", width = 4, height =5 , units = 'in', res = 300)

# create the plot

# Convert distances from centimorgans to megabases
file$dist_MB <- file$dist / 100

# Create a scatter plot with distances in megabases
plot(file$dist_MB, file$rsq, pch=20, col="darkorange", cex=0.5, xlab="Distance (Mb)", ylab="LD (r^2)", cex.axis=1, cex.lab=1.2)

lines(newfile$file.dist, newfile$newrsq, col="navyblue", lwd=4)
abline(h=0.2, col="red", lwd=3, lty="dotted")
abline(v=halfdecaydist, col="green", lwd=3, lty="dotted")
mtext(round(halfdecaydist,2), side=1, line=0.05, at=halfdecaydist, cex=1, col="green")
f1 <- data.frame(newfile$file.dist, newfile$newrsq)
xval <- f1[which.min(abs(0.2 - f1$newfile.newrsq)),] #find x value where y=0.2
xval[1,1]
dev.off()





