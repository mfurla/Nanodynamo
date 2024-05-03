### Libraries
library("ggplot2")

# Untreated - Figure 2A
Unt <- read.table("/path/to/Files/polysomalProfiling/Untreated.csv", header = T, sep = "\t")
df <- data.frame("Reps"=c(rep(1,length(Unt[,"Time1"]))
                          ,rep(2,length(Unt[,"Time2"]))
                          ,rep(3,length(Unt[,"Time3"]))),
                  "Time" = c(Unt[,1],Unt[,3],Unt[,5]),
                  "Absorbance"=c(Unt[,2],Unt[,4],Unt[,6]))
df <- df[-which(is.na(df$Time)),]
df$Time <- as.numeric(df$Time)
df <- df[-which(df$Time>=10.5),]
df <- df[-which(df$Time<=(1.5)),]
# Plot
ggplot(data = df,
             aes(x = Time, y = Absorbance, group=Reps)) + 
  geom_line(aes(color=Reps)) +
  theme_classic() +
  scale_x_discrete(limits = c(1:10)) +  
  xlab("Fraction") +                      
  theme(axis.title.x = element_text(vjust = 0.25),text = element_text(size = 15))  

# Pladienolide B - Figure 3D
PlaB <- read.table("/path/to/Files/polysomalProfiling/PladienolideB.csv", header = T, sep = "\t")
df <- data.frame("Reps"=c(rep(1,length(PlaB[,"Time1"]))
                          ,rep(2,length(PlaB[,"Time2"]))),
                  "Time" = c(PlaB[,1],PlaB[,3]),
                  "Absorbance"=c(PlaB[,2],PlaB[,4]))
df <- df[-which(is.na(df$Time)),]
df$Time <- as.numeric(df$Time)
df <- df[-which(df$Time>=10.5),]
df <- df[-which(df$Time<=(1.5)),]
# Plot

ggplot(data = df,
             aes(x = Time, y = Absorbance, group=Reps)) + 
  geom_line(aes(color=Reps)) +
  theme_classic() +
  scale_x_discrete(limits = c(1:10)) +  
  xlab("Fraction") +                      
  theme(axis.title.x = element_text(vjust = 0.25),text = element_text(size = 15))  

# Leptomycin B - Figure 4A
Lep <- read.table("/path/to/Files/polysomalProfiling/LeptomycinB.csv", header = T, sep = "\t")
df <- data.frame("Reps"=c(rep(1,length(Lep[,"Time1"]))
                          ,rep(2,length(Lep[,"Time2"]))
                          ,rep(3,length(Lep[,"Time3"]))
                          ,rep(4,length(Lep[,"Time4"]))),
                  "Time" = c(Lep[,1],Lep[,3],Lep[,5],Lep[,7]),
                  "Absorbance"=c(Lep[,2],Lep[,4],Lep[,6],Lep[,8]))
df <- df[-which(is.na(df$Time)),]
df$Time <- as.numeric(df$Time)
df <- df[-which(df$Time>=10.5),]
df <- df[-which(df$Time<=(1.5)),]

# Plot
ggplot(data = df,
             aes(x = Time, y = Absorbance, group=Reps)) + 
  geom_line(aes(color=Reps)) +
  theme_classic() +
  scale_x_discrete(limits=c(1:10))+
  xlab("Fraction") +                      
  theme(axis.title.x = element_text(vjust = 0.25),text = element_text(size = 15))  

# Harringtonine - Figure 5A
Harr <- read.table("/path/to/Files/polysomalProfiling/Harringtonine.csv", header = T, sep = "\t")
df <- data.frame("Reps"=c(rep(1,length(Harr[,"Time1"]))
                          ,rep(2,length(Harr[,"Time2"]))),
                  "Time" = c(Harr[,1],Harr[,3]),
                  "Absorbance"=c(Harr[,2],Harr[,4]))
df <- df[-which(is.na(df$Time)),]
df$Time <- as.numeric(df$Time)
df <- df[-which(df$Time>=10.5),]
df <- df[-which(df$Time<=(1.5)),]

# Plot
ggplot(data = df,
             aes(x = Time, y = Absorbance, group=Reps)) + 
  geom_line(aes(color=Reps)) +
  theme_classic() +
  scale_x_discrete(limits = c(1:10)) +  
  xlab("Fraction") +                      
  theme(axis.title.x = element_text(vjust = 0.25),text = element_text(size = 15))