# RSPCA

Code for "Robust Stochastic Principal Component Analysis" by John Goes, Teng Zhang, Raman Arora Gilad Lerman

Paper link: http://proceedings.mlr.press/v33/goes14.pdf

How to use this code:

The following is a table of the estimators used in the paper, both those introduced and those compared via experiments, with names defined in the paper:

SGD - first()

R-SGD1 - firstR1()

R-SGD2 - firstR2()

Inc - second()

R-Inc - secondR()

Mirror Descent (MD) - third()

Robust-Mirror-Descent - thirdR()


The following lists the figures in the paper and which code can reproduce them. Note that typing help function_name() explains the details of how to change parameters, etc.

Figures 1-3: 

compare_all_estimators()  

***Note: There appears to be an error with Figure 1 in the paper. It seems that Figure 1 in the paper corresponds to the setting where 90% of the data are outliers, the inliers live on a needle, and the MD and R-MD algorithms normalized the data to the sphere before estimation.***


Figure 4:

heatmap_gen()

Figure 5:

astrodata_final()

Figure 6:

faces_experiment()

