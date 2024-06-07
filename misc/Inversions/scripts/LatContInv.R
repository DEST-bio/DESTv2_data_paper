args <- commandArgs(TRUE)
library(ggplot2)

##test
#print(args[1])

# Read the CSV file
nhm_inversion <- read.delim(args[1], row.names=1)

for (col in colnames(nhm_inversion)){
  nhm_inversion[[col]] <- as.numeric(paste(nhm_inversion[[col]]))
}  

meta <- read.csv(args[2], row.names=1)
meta_small <- as.data.frame(cbind(meta$lat,as.character(meta$continent)), row.names <- rownames(meta))
colnames(meta_small) <- c("Latitude", "Continent")
joined_table <- merge(meta_small, nhm_inversion, by = "row.names")
data <- joined_table
data <- na.omit(data)

# Define the color variable
color_var <- data$Continent
x <- as.numeric(paste(data$Latitude))
y <- data[[args[3]]]

# Create the ggplot
p <- ggplot(data, aes(y = y, x = x , color = color_var)) +
  geom_point(alpha=.5) +
  labs(x = "Latitude", y = "InversionFrequency", color = "Continent") +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = y ~ x, alpha=.15, aes(fill=Continent)) +
  guides(fill="none") + 
  ggtitle(paste(args[3])) + theme(plot.title = element_text(hjust=0.5)) 


ggsave(paste("/media/DEST2_NHM/results/PoolSNP_",args[3],".pdf", sep=""), p, width=15, height=10)
