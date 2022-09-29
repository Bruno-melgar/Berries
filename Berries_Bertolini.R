##########################################
# Importing Libraries for analysis
##########################################
setwd("~/Desktop")

##########################################
# Importing Libraries for analysis
##########################################
library("drc")
library(tidyverse)
library(ggplot2)
library(patchwork)
rows= 5

##########################################
# attaching main dataset for ease of work
##########################################
attach(deb)

##########################################
# Making selection of important rows
##########################################
rasp1 <- deb %>% 
  select(`Sample ID`, `Sample concentration`, Reduction) %>% 
  filter(`Sample ID` == "TFVR1")

rasp2 <- deb %>% 
  select(`Sample ID`, `Sample concentration`, Reduction) %>% 
  filter(`Sample ID` == "TFVR2")

black1 <- deb %>% 
  select(`Sample ID`, `Sample concentration`, Reduction) %>% 
  filter(`Sample ID` == "TFBla1")

black2 <- deb %>% 
  select(`Sample ID`, `Sample concentration`, Reduction) %>% 
  filter(`Sample ID` == "TFBla2")

jab1 <- deb %>% 
  select(`Sample ID`, `Sample concentration`, Reduction) %>% 
  filter(`Sample ID` == "TFVJ1")

jab2 <- deb %>% 
  select(`Sample ID`, `Sample concentration`, Reduction) %>% 
  filter(`Sample ID` == "TFVJ2")

acai1 <- deb %>% 
  select(`Sample ID`, `Sample concentration`, Reduction) %>% 
  filter(`Sample ID` == "TFVA1")

acai2 <- deb %>% 
  select(`Sample ID`, `Sample concentration`, Reduction) %>% 
  filter(`Sample ID` == "TFVA2")

blue1 <- deb %>% 
  select(`Sample ID`, `Sample concentration`, Reduction) %>% 
  filter(`Sample ID` == "TFVBlu1")

blue2 <- deb %>% 
  select(`Sample ID`, `Sample concentration`, Reduction) %>% 
  filter(`Sample ID` == "TFVBlu2")

straw1 <- deb %>% 
  select(`Sample ID`, `Sample concentration`, Reduction) %>% 
  filter(`Sample ID` == "TFVS1")

straw2 <- deb %>% 
  select(`Sample ID`, `Sample concentration`, Reduction) %>% 
  filter(`Sample ID` == "TFVS2")

##########################################
# Summary for data ceherence 
##########################################
summary(acai1)
summary(acai2)
summary(black1)
summary(black2)
summary(blue1)
summary(blue2)
summary(jab1)
summary(jab2)
summary(rasp1)
summary(rasp2)
summary(straw1)
summary(straw2)

##########################################
# Creating the models
##########################################
modacai1 <- drm(Reduction~`Sample concentration`,
                data = acai1, fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modacai2 <- drm(Reduction~`Sample concentration`,
                data = acai2, fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modblack1 <- drm(Reduction~`Sample concentration`,
                 data = black1, fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modblack2 <- drm(Reduction~`Sample concentration`,
                 data = black2, fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modblue1 <- drm(Reduction~`Sample concentration`,
                data = blue1, fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modblue2 <- drm(Reduction~`Sample concentration`,
                data = blue2, fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modjab1 <- drm(Reduction~`Sample concentration`,
               data = jab1, fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modjab2 <- drm(Reduction~`Sample concentration`,
               data = jab2, fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modrasp1 <- drm(Reduction~`Sample concentration`,
                data = rasp1, fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modrasp2 <- drm(Reduction~`Sample concentration`,
                data = rasp2, fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modstraw1 <- drm(Reduction~`Sample concentration`,
                 data = straw1, fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modstraw2 <- drm(Reduction~`Sample concentration`,
                 data = straw2, fct = LL.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

##########################################
# Summarizing the models for modeling coherence
#The estimate e: its the EC50 value
##########################################
summary(modacai1)
#summary(modacai2)
#summary(modblack1)
summary(modblack2)
summary(modblue1)
#summary(modblue2)
#summary(modjab1)
summary(modjab2)
summary(modrasp1)
#summary(modrasp2)
summary(modstraw1)
#summary(modstraw2)

#interval = "delta" gives you asymptotically-based confidence intervals at a default 95% level. No funciona con Weibull 
ED(modacai1W2.4, c(10,20,50), interval="delta")

##########################################
# choosign the model
#Lower Residual and Akaike´s information criterion (IC) the better
# Lack of fit, higher the p-value the better
##########################################
model.LL3<- drm(Reduction~`Sample concentration`,
                data = straw1, fct=LL.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit", "ED50")))
mselect(model.LL3, fctList = list(W1.3(fixed=c(NA, 100, NA)),
                                  W1.4(), 
                                  W2.3(fixed=c(NA, 100, NA)), 
                                  W2.4(),  
                                  LL.4()), linreg=TRUE)

#Best models selected
modacai1W2.4 <- drm(Reduction~`Sample concentration`,
                    data = acai1, fct = W2.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modblack2W2.4 <- drm(Reduction~`Sample concentration`,
                     data = black2, fct = W2.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modblue1W2.4 <- drm(Reduction~`Sample concentration`,
                    data = blue1, fct = W2.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modjab2W2.4 <- drm(Reduction~`Sample concentration`,
                   data = jab2, fct = W2.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modrasp1W2.4 <- drm(Reduction~`Sample concentration`,
                    data = rasp1, fct = W2.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

modstraw1W2.4 <- drm(Reduction~`Sample concentration`,
                     data = straw1, fct = W2.4(names = c("Slope", "Lower limit", "Upper limit", "EC50")))

# Plot Comparison of best models
plot(modstraw1, type = "all")
plot(modstraw1W2.4, add = TRUE, type = "all", col= "orange")

#Summary of the best model
a1<- summary(modacai1W2.4)
b1<- summary(modblack2W2.4)
c1<- summary(modblue1W2.4)
d1<- summary(modjab2W2.4)
e1<- summary(modrasp1W2.4)
f1<- summary(modstraw1W2.4)

sink("f.txt")
f1
sink()

##########################################
# checking singular plots
##########################################
plot(modacai1, type = "all", xlab = "[mg/mL]", ylab = "% Scavenging", main= "Açai Maceration")
plot(modblack2, type = "all", xlab = "[mg/mL]", ylab = "% Scavenging", main= "Blackberry Maceration")
plot(modblue1, type = "all", xlab = "[mg/mL]", ylab = "% Scavenging", main= "Blueberry Maceration")
plot(modjab2, type = "all", xlab = "[mg/mL]", ylab = "% Scavenging", main= "Jabuticaba Maceration")
plot(modrasp1, type = "all", xlab = "[mg/mL]", ylab = "% Scavenging", main= "Raspberry Maceration")
plot(modstraw1, type = "all", xlab = "[mg/mL]", ylab = "% Scavenging", main= "Strawberry Maceration")

##########################################
# Creating multiple plot
##########################################
par(mar=c(5,4,4,10), xpd=TRUE)
plot(modacai1W2.4, type = "all", xlab = "Sample concentration [mg/mL]", ylab = "DPPH Scavenging Effect (%)", 
     col = "blueviolet", pch=1, #legend=TRUE, legendText = "Açai", 
     lty = 1, lwd =2, xlim = c(0.01,3.2), ylim = c(0,100))
plot(modblack2W2.4, add = TRUE, type = "all",
     col = "black", pch=2, #legend=TRUE, legendText = "Blackberry", 
     lty = 1, lwd =2)
plot(modblue1W2.4, add = TRUE, type = "all",
     col = "blue", pch=3, #legend=TRUE, legendText = "Blueberry",
     lty = 1, lwd =2)
plot(modjab2W2.4, add = TRUE, type = "all",
     col = "darkgreen", pch=4, #legend=TRUE, legendText = "Jabuticaba",
     lty = 1, lwd =2)
plot(modrasp1W2.4, add = TRUE, type = "all",
     col = "magenta", pch=5, #legend=TRUE, legendText = "Raspberry",
     lty = 1, lwd =2)
plot(modstraw1W2.4, add = TRUE, type = "all",
     col = "red", pch=6, #legend=TRUE, legendText = "Strawberry",
     lty = 1, lwd =2)
legend("topright", inset = c(-0.2,0), 
       legend = c("Açai", "Blackberry", "Blueberry", "Jabuticaba", "Raspberry", "Strawberry"), 
       pch = c(1,2,3,4,5,6), text.col = c("blueviolet","black","blue","darkgreen","magenta","red"),
       title = "Fruits")


##########################################
# attaching main dataset for ease of work and data manipulation
##########################################
attach(antiox)

antiox$Sample <- as.factor(antiox$Sample)

summary(antiox)

##########################################
# Creating Facet plot
##########################################
antiox %>% 
  ggplot(aes(`Antiox assay`,`IC50 (mg/mL)`, fill=`Antiox assay`)) +
  geom_bar(stat='identity',) + 
  facet_wrap(~Sample, ncol = 6) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0))
#theme_minimal()