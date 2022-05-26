

#Author: Annalisa Bellandi, Faulkner group (John Innes centre)

#Title: Linear vascular progression analysis

#Content: given data of linear distance from the wounded site over time
#allows for calculation of the velocity in 2 ways
#1. fitting a polynomial of degree 3 and calculating its derivative
#2. fitting a stright line and extracting the slope

#In the case of progression of calcium waves along the vasculature, the disposition of the measured datapoints suggests that 
#a straight line is a better fit for the phenomenon.
#Slopes of the lines are therefore plotted in boxplots
#comaprison across series (genotypes/conditions) is achieved by wilcoxon rank sum test


##-------------------------------------------------------------------------------------------------------------------------
##------------------------------------------ packages needed  --------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------

library("ggplot2")
library('ggsignif')
library('ggpubr')
library('rlang')
library("purrr")
library("tidyr")
library("Rmisc")
library("reshape2")
library("dplyr")
library("openxlsx")
library("RColorBrewer")
library("svglite")

##-------------------------------------------------------------------------------------------------------------------------
##------------------------------------------ set my colour palette --------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------

my_pal_div <- RColorBrewer::brewer.pal(11, "BrBG")[2:11]
my_pal_quant_1 <- RColorBrewer::brewer.pal(9, "Oranges")
my_pal_quant_2 <- RColorBrewer::brewer.pal(9, "Blues")
my_pal_gray <- RColorBrewer::brewer.pal(9, "Greys")
okabe_pal <- c("#E69F00","#56B4E9","#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7")

n <- max(length(my_pal_div), length(my_pal_quant_1), length(my_pal_quant_2), length(my_pal_gray), length(okabe_pal))

length(my_pal_div) <- n
length(my_pal_quant_1) <- n
length(my_pal_quant_2) <- n
length(my_pal_gray) <- n
length(okabe_pal) <- n

my_pal_gray_d <- data.frame(my_pal_gray)
my_pal_quant_1_d <- data.frame(my_pal_quant_1)
my_pal_quant_2_d <- data.frame(my_pal_quant_2)
my_pal_div_d <- data.frame(my_pal_div)
okabe_pal_d <- data.frame(okabe_pal)

my_col_set <- (0)
my_col_set <- cbind(my_pal_gray_d, my_pal_quant_1_d)
my_col_set <- cbind(my_col_set, my_pal_quant_2_d)
my_col_set <- cbind(my_col_set, okabe_pal_d)
my_col_set <- cbind(my_col_set, my_pal_div_d)

my_col_set_df <- data.frame(my_col_set)

order <- c(1:10)
my_col_set_df1 <- cbind(my_col_set_df, order)
my_col_set_df1

long_color <- melt(my_col_set_df1,
                   id.vars = "order",
                   variable.name = "palette",
                   value.name = "color")

my_colors_plot <- ggplot(long_color, aes(x = palette, y = order, fill = color)) +
  geom_tile(aes(width=0.93, height=0.95)) +
  scale_fill_identity() +
  scale_y_continuous(breaks=c(1:n)) +
  theme_light()+
  geom_label(aes(label=color),colour = "black", fill= "white", fontface = "bold", size = 4)

my_colors_plot


##-------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------


##=======================================================================================================================
##=======================================  to manually set before start  ====================================================

setwd ("...") #set working dyrectory
GEN='...' #set the name the data series that you will use (genotype)
CON='...' #set another layer of information on the series if needed (eg condition)

#From the excel file, copy data-frame with data of progression of the front of the signal over time, for one genotype/condition across all replicates
df <- read.table('clipboard', sep = '\t', header=TRUE, stringsAsFactors = FALSE)
head(df,10)
df<-data.frame(df)
tail(df)


##=========================================== analysis of the velocities of progression of each replicate for this genotype/condition ============================================
##===================================== this chunk need to be re-run separately for each genotype/condition of one experiment 

rep <- unique(df$series)

rep

fit_poly <- data.frame(genotype=character(), sample_name=numeric(), condition=character(), coeff1=numeric(), coeff2=numeric(), coeff3=numeric())

slopes <- data.frame(genotype=character(), sample_name=numeric(), condition=character(), alpha=numeric())

workbook_fits <- createWorkbook(paste(GEN, 'fits', ".xlsx"))


#for each replicate retrieves the postion fo thr front of the wave over time 

for (j in rep) {
  
  df_rep <- df[df$series==j,]
  
  #---------------plot the data in scatter plot
  
  base_plot <- ggplot()+
    geom_point(data=df_rep, aes(x=time, y=d), shape=21, fill=okabe_pal[1], col=my_pal_gray[5], stroke=0.01, size=3.5, alpha=1)+
    xlab("Time [s]") +
    ylab(paste("Distance from the wound", "\n" , "[\U003BCm]"))+
    theme_bw(base_size=18)+
    theme(plot.title = element_text(hjust = 0.5))+
    ggtitle(paste("Signal along the main vein", '\n', GEN, "series", j, '\n', CON, sep=""))
  
  
  base_plot
  
  #saves PNG version for easy acces
  ggsave(paste(GEN, "series", j, "distance vein", ".png"), base_plot, width = 8, height = 5, dpi = 450, device = 'png')
  
  
  
  #-----------------fits a polynomial of degree 3 on the data and plots it
  
  model_poly3 <- lm(df_rep$d ~ poly(df_rep$time,3, raw = T) -1) # the -1 forces the poly to pass from 0,0
  cf0 = coef(model_poly3)
  cf0
  
  xm <- 0:max(df_rep$t)
  y0 <- cf0[1]*xm + cf0[2]*(xm^2) + cf0[3]*(xm^3) #+ cf0[4]*(xm^4) #+ cf0[5]*(xm^5) + cf0[6]*(xm^6) #could increase the degrees of the polynomial if needed
  model_poly6_df <- data.frame(cbind(xm,y0))
  
  
  coefficient_poly <- data.frame(genotype=GEN, sample_name=j, condition=CON)
  
  coefficient_poly$coeff1 <- cf0[1]
  coefficient_poly$coeff2 <- cf0[2]
  coefficient_poly$coeff3 <- cf0[3]
  #coefficient_poly$coeff4 <- cf0[4]
  #coefficient_poly$coeff5 <- cf0[5]
  #coefficient_poly$coeff6 <- cf0[6]
  
  coefficient_poly
  
  coefficient_poly3 <- coefficient_poly
  
  #now plot the model line on top and then save
  
  poly3_plot <- base_plot +
    geom_line(data=model_poly6_df, mapping=aes(x=model_poly6_df$xm, y=model_poly6_df$y0),
              color=my_pal_gray[9], lty=1, lwd=1.5, alpha=0.3)
  
  poly3_plot
  
  #saves PNG version for easy acces
  ggsave(paste(GEN, "series", j, "polynomial constrained fit 3rd", ".png"), poly3_plot, width = 8, height = 5, dpi = 450, device = 'png')
  
  
  
  #----------------- first derivative of the polynomial model - this gives us the velocity at eavry time point
  xmD1 <- 0:max(df_rep$t)
  
  D1 = cf0[1] + 2*cf0[2]*xmD1 + 3*cf0[3]*xmD1^2 #+ 4*cf0[4]*xmD1^3 + 5*cf0[5]*xmD1^4 + 6*cf0[6]*xmD1^5
  
  model3tD1_df <- data.frame(cbind(xmD1,D1))
  
  poly3_D1 <- ggplot()+
    geom_line(data=model3tD1_df, mapping=aes(x=xmD1, y=D1),
              color=my_pal_gray[9], lty=1, lwd=1.5, alpha=0.3) +
    xlab("Time [s]") + 
    ylab("Velocity [\U003BCm/s]") + 
    ggtitle (paste("Velocity along main vein", '\n', GEN, "series", j, '\n', CON)) + 
    theme_bw(base_size=18)+
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.position="none")
  
  
  poly3_D1
  
  #saves PNG version for easy acces
  ggsave(paste(GEN, "series", j, "derivative poly 3th", ".png"), poly3_D1, width = 8, height = 5, dpi = 450, device = 'png')
  
  fit_poly <- rbind(fit, coefficient_poly3)
  
  
  #----------------- fits a straight line on the data and plots it
  
  model_line <- lm(df_rep$d ~ df_rep$time -1)
  coeff_model = coef(model_line)
  coeff_model
  
  slope_df <- data.frame(genotype=GEN, sample_name=j, condition=CON)
  
  slope_df$alpha <- coeff_model[1]
  
  slopes <- data.frame(rbind(slopes, slope_df))
  
  y0=coeff_model[1]*xm
  model_line_df <- data.frame(cbind(xm,y0))
  
  
  line_plot <- base_plot +
    geom_line(data=model_line_df, mapping=aes(x=xm, y=y0),
              color=my_pal_gray[9], lty=1, lwd=1.5, alpha=0.3)
  
  
  line_plot
  
  #saves PNG version for easy acces
  ggsave(paste(GEN, "series", j, "line fit", ".png"), line_plot, width = 8, height = 5, dpi = 450, device = 'png')
  
}


addWorksheet(workbook_fits, '3rd_poly')
writeData(workbook_fits, '3rd_poly', fit_poly)
saveWorkbook(workbook_fits, file=paste(paste(GEN, 'fits', ".xlsx")), overwrite = TRUE)

addWorksheet(workbook_fits, 'line')
writeData(workbook_fits, 'line', slopes)
saveWorkbook(workbook_fits, file=paste(paste(GEN, 'fits', ".xlsx")), overwrite = TRUE)


xm <- 0:30
df_reps <- c()
steps <- rownames(slopes)

for (i in steps) {
  y <- slopes[i,4]*xm
  replicate <- cbind(xm, y)
  df_reps <- rbind(df_reps, replicate)
}

head(df_reps)

df_reps <- data.frame(df_reps)

df_reps <- cbind(df_reps, GEN)

head(df_reps,10)

tail(df_reps,10)

#summarySE requires Rmisc
summ_df <- summarySE(data=df_reps, measurevar="y", groupvars="xm")

head(summ_df,10)

summ_df <- cbind(summ_df, GEN)

head(summ_df,10)


#---------------- to show the variability between replicates, plot all the straight fit lines together

plot_all <-  ggplot()+
  geom_point(data=df_reps, mapping=aes(x=xm, y=y),
             color=my_pal_gray[9], stroke=0.01, size=3.5, alpha=0.3) +
  xlab("Time [s]") +
  ylab(paste("Distance from the wound", "\n" , "[\U003BCm]"))+
  theme_bw(base_size=18)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle(paste("Fitted lines,", " all replicates", sep=""))

plot_all

#saves PNG version for easy acces
ggsave(paste(GEN, "line fit_all", ".png"), plot_all, width = 8, height = 5, dpi = 450, device = 'png')

ggsave(paste(GEN, "line fit_all", ".pdf"), plot_all, width = 8, height = 5, dpi = 450)


addWorksheet(workbook_fits, 'line_reps')
writeData(workbook_fits, 'line_reps', df_reps)
saveWorkbook(workbook_fits, file=paste(paste(GEN, 'fits', ".xlsx")), overwrite = TRUE)



#---------------- plot the average line across all replicates

avg_plot <- ggplot()+
  geom_line(data=summ_df, mapping=aes(x=xm, y=y),
            color=my_pal_gray[9], lty=1, lwd=1.5, alpha=0.3) +
  geom_ribbon(data=summ_df, mapping=aes(x=xm, ymin=y+se, ymax=y-se), alpha=0.1) +
  xlab("Time [s]") +
  ylab(paste("Distance from the wound", "\n" , "[\U003BCm]"))+
  theme_bw(base_size=18)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle(paste("Average across all replicates", sep=""))

ggsave(paste(GEN, "average_line", ".png"), avg_plot, width = 8, height = 5, dpi = 450, device = 'png')

ggsave(paste(GEN, "average_line", ".pdf"), avg_plot, width = 8, height = 5, dpi = 450)


addWorksheet(workbook_fits, 'line_summ')
writeData(workbook_fits, 'line_summ', summ_df)
saveWorkbook(workbook_fits, file=paste(paste(GEN, 'fits', ".xlsx")), overwrite = TRUE)





##==============================================================================================================================================
#====================================== RUN until here for each genotype/condition separately ==================================================

#============================= then copy all the data frames of different genotypes in one sheet and then run the following for total plots ===================


#------------------------------------------------------ to set before start ----------------------------
#-------------------------------------------------------------------------------------------------------

setwd ("...") #set working dyrectory
EXP='...' #set the name the experiment

GENO1 <- '...' #insert genotype 1 as it is named in the datasets
GENO2 <- '...' #insert genotype 2 as it is named in the datasets

#copy the total summ_df - data-frame with data of progression of the front of the signal over time
df <- read.table('clipboard', sep = '\t', header=TRUE, stringsAsFactors = FALSE)
head(df,10)
df<-data.frame(df)
tail(df)

#copy the total data frame that contains the slopes of each replicate of all genotypes/sereis/condition
df_slopes <- read.table('clipboard', sep = '\t', header=TRUE, stringsAsFactors = FALSE)
head(df_slopes,10)
df_slopes<-data.frame(df_slopes)
tail(df_slopes)



#------------------------------------------this plots the average line for each genotype/series/condition 

avg_plot <- ggplot(data=df, mapping=aes(x=xm, y=y, col=GEN))+
  geom_line(lty=1, lwd=1.5, alpha=0.8)+
  geom_ribbon(aes(x=xm, ymin=y+se, ymax=y-se, fill=GEN, col=NA), alpha=0.2) +
  scale_color_manual(values=c(okabe_pal[1:7]))+
  scale_fill_manual(values=c(okabe_pal[1:7]))+
  xlab("Time [s]") +
  ylab(paste("Distance from the wound", "\n" , "[\U003BCm]"))+
  theme_bw(base_size=18)+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle(paste("Average across all replicates", sep=""))+
  theme(text =element_text(size=20)) +
  theme(#legend.key = element_rect(fill = NA),
    legend.key.size = unit(1.5, "cm"), legend.position = 'bottom')

avg_plot

ggsave(paste(EXP, "average_line_all", ".png"), avg_plot, width = 8, height = 6, dpi = 450, device = 'png')

ggsave(paste(EXP, "average_line_all", ".pdf"), avg_plot, width = 8, height = 6, dpi = 450)





#------------------------------------------ this plots boxplots for the slopes with wilcoxon rank sum test


plot_wilcox <- ggplot()+
  geom_boxplot(df, mapping=aes(x=genotype, y=alpha, color=genotype), outlier.shape = NA)+
  geom_point(df, mapping=aes(x=genotype, y=alpha, color=genotype, fill=genotype), 
             size = 3, shape = 21, alpha=0.4, position = position_jitterdodge()) +
  scale_color_manual(values = c(okabe_pal[1], okabe_pal[2]))+
  scale_fill_manual(values = c(okabe_pal[1], okabe_pal[2]))+ #dd my_pal_gray[5] if needed col0
  xlab("genotype")+
  theme_bw(base_size=28)+ 
  #coord_cartesian(ylim=c(0,500))+
  guides(shape=FALSE)+
  geom_jitter(width=0.1, alpha=.4)+
  ylab(paste("Velocity [\U003BCm/s]"))+
  ggtitle("Velocity along main vein")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")+
  stat_summary(df, mapping = aes(x=genotype, y=alpha, color=genotype), fun.y=mean, geom="point", shape=4, size=4, stroke=2)+
  geom_signif(data=df, comparisons =list(c(GENO1, GENO2)),
              mapping = aes(x=genotype, y=alpha),
              step_increase = 0.1, map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05), test = 'wilcox.test', color='black',size = 1, y_position = 450, textsize = 7) #default is to wilcoxon test



#saves PNG version for easy acces
ggsave(paste(EXP, "slopes_wilcox", ".png"), plot_wilcox, width = 8, height = 8, dpi = 450, device = 'png')

#saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
#strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
svglite(file=paste(EXP, "slopes_wilcox", ".svg"),width=8,height=8)
figure <- plot_wilcox
print(figure)
dev.off()

#saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
#illustrator, however this renders correctly in inkscape
ggsave(paste(EXP,"slopes_wilcox", ".pdf"), plot_wilcox, width = 8, height = 8, dpi = 450)




#------------------------------------------ this calculates separate statistics, in case the exact p value is needed

#requires ggpubr
#calculates p value of comaprison based on wilcoxon rank sum test
compare_means(alpha ~ genotype, data = df_slopes, paired = F)

