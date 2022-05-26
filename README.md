# Linear-progression

This repository contains the R script to analyze the progression of a signal along a linear path over time.

Starting data are obtained by image analysis in Fiji, using the straight line tool to measure the distance from the starting site reached by the signal at successive time points. Information on the distance and on the time of each measure can be saved in an excel file with columns:
 - series: the number of that particular sample
 - time: the time of the measurament, relative to teh start of the signal
 - distance: the distance measured by the wave in that time

Create one excel file containing all samples for each genotype/condition. These excel files can then be used for the analysis in R following the script provided here.
