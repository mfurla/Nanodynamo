# Code for Figure S6

# For each txt file in Files/polysomalProfiling

df <- read.table("/path/to/txt file", header = T, sep = "\t")
df <- df %>%  
  select(c(1,2)) %>%      #keep only time and absorbance
  rename(Time = Time..min... , AU = UV..A.U...) %>%   #rename columns
  filter(Time >= 1.500 & Time <= 10.5) #descard non infromative values 

# Plot
g1 <- ggplot(data = df,
             aes(x = Time, y = AU)) + 
  geom_line() +
  theme_classic() +
  scale_x_discrete(limits = c(1:10)) +  
  xlab("Fraction") +                      
  theme(axis.title.x = element_text(vjust = 0.25),text = element_text(size = 15))  

