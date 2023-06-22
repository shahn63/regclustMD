# regclustMD
################################################################################
##
##   README accompanying the publication
##
##   Model-based clustering of mixed data with sparse dependence
##
##   by Young-Gen Choi<ygchoi@skku.edu>, 
##       Soohyun Ahn<shahn@ajou.ac.kr>, and 
##       Jayoun Kim<nunadli03@snu.ac.kr>.
##
################################################################################

########################################
# Description of content:
########################################

* ar1_model.R:
  -------

  Code to reproduce the results of AR1 model in Table 2 presented in the manuscript.

  - One needs to load the main source code 'solver1.54.R' e.g. via

      source("solver1.54.R")
      sourceCpp("cal_prob.cpp") 

  The code with several parameters produces table 2 of simulation study in Section 6.   

* random_model.R:
  -------

  Code to reproduce the results of Random model in Table 2 presented in the manuscript.

  - One needs to load the main source code 'solver1.54.R' and 'cal_prob.cpp' e.g. via

      source("solver1.54.R")
      sourceCpp("cal_prob.cpp") 

  The code with several parameters produces table 2 of simulation study in Section 6. 

* Prostate.R:
  
  Code to reproduce the analysis of the Prostate cancer data in Table 7 presented in the manuscript.

  - One needs to load the main source code 'solver1.54.R', 'cal_prob.cpp', and 'cluster_validation_final.R' e.g. via

      source("solver1.54.R")
      sourceCpp("cal_prob.cpp") 
      source("cluster_validation_final.R") 

* AIS.R:
  
  Code to reproduce the analysis of the Australian Institute of Sports data in Table 8 presented in the manuscript.

  - One needs to load the main source code 'solver1.54.R', 'cal_prob.cpp', and 'cluster_validation_final.R' e.g. via

      source("solver1.54.R")
      sourceCpp("cal_prob.cpp") 
      source("cluster_validation_final.R") 

* README.txt:
  -----------

  This file.

########################################
# R programming specifications
########################################

Session info:

It was tested with the following configuration:

    R version 4.0.3 (2020-10-10)
    Platform: x86_64-pc-linux-gnu (64-bit)


