##### Script to visualize lab results ########
## Clear directory ##
rm(list = ls())

## Import necessary libraries ##
library(tidyverse)
library(reshape2)
library(plyr)


# ### Import data ###
# cfu_data <- read.csv(file = "C:/Users/Zachary/Documents/Edinburgh/Lab_results/cfu_counts_20190717_Pento.csv",
#                      header = T)
# 
# head(cfu_data)
# str(cfu_data)
# 
# plot_data <- cfu_data[,c(1,2,4,6,8)]
# avg_plot_data <- ddply(plot_data, .(treatment), summarize, mean_cfus_day1 = mean(day1_log_cfus, na.rm = T),
#       mean_cfus_day2 = mean(day2_log_cfus, na.rm = T), mean_cfus_day3 = mean(day3_log_cfus, na.rm = T))
# avg_plot_data <- avg_plot_data[c(-5,-6),]
# 
# avg_plot_data <- melt(data =avg_plot_data, value.name = 'log_cfus')
# 
# head(avg_plot_data)
# 
# ggplot(data = avg_plot_data, aes(x = variable, y = log_cfus)) +
#   geom_jitter(aes(color = treatment), size = 2.4) +
#   labs(x = "Days after oral infection", y= "Log10 CFUs")
# 


#### Import data from last experiment ####

cfu_data <- read.csv(file = "C:/Users/Zachary/Documents/Edinburgh/Lab_results/cfu_counts_pento_20190802_complete.csv",
                     header = T)
head(cfu_data)
str(cfu_data)

plot_data <- ddply(cfu_data, .(treatment, day), summarize, mean_cfus = mean(cfus, na.rm = T),
                   mean_log10_cfus = mean(log10_cfus, na.rm = T),
                   se_log10_cfus = (sd(log10_cfus, na.rm = T)/sqrt(length(log10_cfus[!is.na(log10_cfus)]))))
colnames(plot_data)[1] <- "Treatment"

### Plot of the average number of Log CFUs in each treatment with zeros included ### 
ggplot(data = plot_data, aes(x = day, y = mean_log10_cfus)) +
  geom_line(aes(color = Treatment), position = position_dodge(0.08))+
  geom_point(aes(color = Treatment), size = 2, position = position_dodge(0.08)) +
  labs(x = "Days after Oral Infection", y = "Mean Log10 CFUs by Treatment") +
  geom_errorbar(aes(ymin = mean_log10_cfus - se_log10_cfus, ymax = mean_log10_cfus + se_log10_cfus, color = Treatment), 
                width = 0.1, position = position_dodge(0.08))

### Removing zeros ###
cfu_data_zero <- cfu_data[which(cfu_data$log10_cfus!=0),]

cfu_data[which(cfu_data$log10_cfus!=0),]
length(which(cfu_data$log10_cfus==0 & cfu_data$day == 1))
length(which(cfu_data$log10_cfus==0 & cfu_data$day == 2))
length(which(cfu_data$log10_cfus==0 & cfu_data$day == 3))
length(which(cfu_data$log10_cfus==0 & cfu_data$day == 4))
length(which(cfu_data$log10_cfus==0 & cfu_data$day == 5))


plot_data <- ddply(cfu_data_zero, .(treatment, day), summarize, mean_cfus = mean(cfus, na.rm = T),
                   mean_log10_cfus = mean(log10_cfus, na.rm = T),
                   se_log10_cfus = (sd(log10_cfus, na.rm = T)/sqrt(length(log10_cfus[!is.na(log10_cfus)]))))

### Plot of the average number of Log CFUs in each treatment without zeros included ### 
ggplot(data = plot_data, aes(x = day, y = mean_log10_cfus)) +
  geom_line(aes(color = treatment))+
  geom_point(aes(color = treatment), size = 2, position = position_dodge(0.05)) +
  labs(x = "Days after Oral Infection", y = "Mean Log10 CFUs by Treatment") +
  geom_errorbar(aes(ymin = mean_log10_cfus - se_log10_cfus, ymax = mean_log10_cfus + se_log10_cfus, color = treatment), 
                width = 0.1, position = position_dodge(0.05))

#### A statistical analysis of our CFU counts ####




#### Comparing the MICs of our P ento samples #####

# Import data #
df_mic <- read.csv(file = "C:/Users/Zachary/Documents/Edinburgh/Lab_results/mic_plate_results_pento.csv",
         header = T, strip.white = T)

head(df_mic)

str(df_mic)

df_mic <-df_mic[which(is.na(df_mic$mic)==FALSE),]

df_mic_sub <- subset(df_mic, treatment != "1% Methanol" & treatment != "1.5% Methanol")
df_mic_sub$treatment <- droplevels(df_mic_sub$treatment)

df_mic_totals <- ddply(df_mic_sub, .(treatment, resistance), .fun = summarize, resistance_total = length(sample_id))

df_mic_wide <- spread(df_mic_totals, key = treatment, value = resistance_total)
table(df_mic_wide[1,])

#### Plotting the MIC totals of our flies ####
df_mic_totals$resistance <- factor(x = df_mic_totals$resistance, levels = c("susceptible","resistant"), ordered = TRUE)

colnames(df_mic_totals)[1] <- "Treatment"
ggplot(df_mic_totals, aes(x = resistance, y= resistance_total)) + 
  geom_bar(stat = "identity", aes(fill = Treatment)) +
  labs(x = "Antimicrobial Resistance Status", y= "Number of Fly Derived Bacteria Samples")
  
as.matrix(df_mic_wide)
df_mic_wide[2,2] <- 0
df_mic_wide[2,3] <- 0

### Split this into two comparisons the 100 and 150 ug/ml treatments versus the untreated
matrix_100 <- as.matrix(df_mic_wide[,c(2,4)])

fisher.test(x = matrix_100, alternative = "greater")

matrix_150 <- as.matrix(df_mic_wide[,c(3,4)])

fisher.test(x = matrix_150, alternative = "greater")


prop.test(df_mic_totals$mic_total[11:18], grp_sz,  alternative = "greater")

chisq.test(df_mic$treatment[c(1:10,21:30)], df_mic$mic[c(1:10,21:30)])
table(df_mic$treatment, df_mic$mic)
