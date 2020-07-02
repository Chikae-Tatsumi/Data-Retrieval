# Made by Yoshitaka Uchida & Chikae Tatsumi

setwd("~/R/Analysis/1_Test")

# Import and enter your data!!!!
data <- read.csv("admin_XXXXX -  Quantification Amplification Results_SYBR.csv", header = TRUE) # Replace fine name
soil <- read.csv ("soil data.csv",header=T, row.names=1) # "weight_g" column = extracted fresh soil weight (g), "water content_%" column = water content of each soil (%), "DNA_ngperml" column = DNA quantity from the quanttus or something (ng/ml) # see example file
std.conc <- 8144  # ng/ul, from Qubit or something like that
std.dilt <- c(10^-4,10^-5,10^-6,10^-7,10^-8,10^-9)   # How did you dilute your standard samples? 
std.number <- 6 # How many wells did you use for the standard samples?
nc.number <- 2 # How many wells did you use for the negative control?
std.rm <- c(3,6)  # or c() # If you want to remove some of the standard samples from the calculation, what number of the standand samples do you want to remove?
sample.dilt <- 50 # How much did you dilute the samples?
DNA.length <- 291 # bp, How long is your stndrad DNA? if you use the amplified 16S genes, it's 291 bp. See http://omegabioservices.com/index.php/16s-reference/
final.nfw <- 30 # ul, how much of nuclease-free water was used to extract DNA at the final step


# Get ready
if(!require("svDialogs")){install.packages("svDialogs")}
library(svDialogs)
if(!require("qpcR")){install.packages("qpcR")}
library(qpcR)
if(!require("stringr")){install.packages("stringr")}
library(stringr)
if(!require("ggplot2")){install.packages("ggplot2")}
library(ggplot2)


# Basic calculation
DNA.weight.mol <- DNA.length * 330 *2  # g/mol
DNA.weight <- DNA.weight.mol/(6.023*10^23)  # g/molecule
DNA.amount <- 10^(-6)/DNA.weight  # molecules/ug


# Illustrate the amplification curves and get the Ct values
# data <- data[,-1] # Delete the first column
data <- data[,c(1,2,14,26,38,50,62,74,86,3,15,27,39,51,63,75,87,4,16,28,40,52,64,76,88,5,17,29,41,53,65,77,89,6,18,30,42,54,66,78,90,7,19,31,43,55,67,79,91,8,20,32,44,56,68,80,92,9,21,33,45,57,69,81,93,10,22,34,46,58,70,82,94,11,23,35,47,59,71,83,95,12,24,36,48,60,72,84,96,13,25,37,49,61,73,85,97)]
m1 <- modlist(data, 1, 2:97, l4)

png(file = "results.png") # Save the amplification curves
par(ps = 18)
plot(m1, type = c("single"))
dev.off()

m2 <- getPar(m1, type = c("curve"))
m2 <- t(m2)
head(m2)


# Make the standard curve
m2 <- data.frame(m2)
std.data <- m2[1:(std.number),]
std.data$conc <- std.dilt*std.conc 
std.data$copy.number <- (std.data$conc*DNA.amount)
std.data$log.copy.number <- log10(std.data$copy.number)
std.data <- std.data[-(std.rm),]

std.curve <- lm(ct ~ log.copy.number, data = std.data)  
intercept <- summary(std.curve)[[4]][1,1] 
slope <- summary(std.curve)[[4]][2,1]
r2.value <- summary(std.curve)[[9]] # Adjusted R-squared
if (r2.value < 0.98) 
{print("Your standard curve is not so good...You should reconsider the standard samples to be used for the calculation of the standard curve")}else{print("Your standard curve is perfect!!! Conglatulations!!!")}
ggplot (std.data, aes(x=log.copy.number,y=ct))+
geom_point()+
geom_smooth(formula = y ~ x,method = "lm",col="black") +
annotate("text", label= paste("slope = ", slope),x=Inf,y=-Inf,hjust=5.5,vjust=-7)+
annotate("text", label= paste("intercept = ", intercept),x=Inf,y=-Inf,hjust=2.3,vjust=-5)+
annotate("text", label= paste("R2 = ", r2.value),x=Inf,y=-Inf,hjust=2.5,vjust=-3)
ggsave ("standard curve.png")


# Calculation of copy numbers
result <- m2[(std.number+nc.number+1):(std.number+nc.number+nrow(soil)),] # get sample wells
result <-cbind (result, soil)

result$conc <- (10^((result$ct-intercept)/slope))*sample.dilt  # copies/ul

result$conc.extraction <- result$conc*final.nfw  # copies in extraction
result$soil.dry.weight <- result$weight_gram*(1-(result$water_content_percent/100)) #g soil
result$copies.gsoil <- result$conc.extraction/result$soil.dry.weight # copies/g soil
result$log.copies.gsoil <- log10(result$copies.gsoil) # log

result$copies.ngDNA <- result$conc/(result$DNA_ngperml/1000) # copies/ng DNA
result$log.copies.ngDNA <- log10(result$copies.ngDNA) # log


# Save
write.csv(result, "qPCR result.csv")