##The version of R needs to be 4.1.2 or later.

##Rtools needs to be downloaded and installed in advance.
##You can download and install the latest Rtools at the official website (https://cran.r-project.org/bin/windows/Rtools/).

##Install the latest version of "cmdstanr" package.
##The "cmdstanr" package can not be installed by CRAN, and the installation process of "cmdstanr" is as following:
##Step1: install "git" software which can download at the website (https://git-scm.com/)
##Step2: set the path of Rtools "C:\rtools40\mingw64\bin", Rtools "C:\rtools40\usr\bin" and git "C:\Program Files\Git\cmd" into "advanced System Settings - environment variable - system variable - path" of your computer.
##Step3: install the cmdstanr" package
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
##Step4: install the "cmdstan" working environment
##Check your toolchain is set up properly
library(cmdstanr)
check_cmdstan_toolchain()
##The correct answer is "The C++ toolchain required for CmdStan is setup properly!". If error happens, see website (https://mc-stan.org/docs/cmdstan-guide/cmdstan-installation.html)
##install cmdstan
install_cmdstan()
##To check the path to the CmdStan installation and the CmdStan version number you can use cmdstan_path() and cmdstan_version():
cmdstan_path()
cmdstan_version()

##Install the latest version of "kinship2" package.
install.packages("kinship2")

##Install the "MBEXCIS" package.
##setwd() to the root of the folder path of "BEXCIS_1.0.zip" and run the following command
install.packages("MBEXCIS_1.0.zip",repos = NULL)
##or setwd() to the root of the folder path of "BEXCIS_1.0.tar.gz" and run the following command
install.packages("MBEXCIS_1.0.tar.gz",type='source')

###########################################################################
##################    examples for M_Bayes_XCI()    #########################
###########################################################################

##example 1:
##There is a data set "mixture_data1" that includes a quantitative trait and the genotype of target SNP of 30 pedigrees, 260 pedigree-related individuals and extra 130 unrelated females.
##There is the corresponding covariate data "covariate_mixture1_2" for above individuals.

############################  step 1: loading the R package  ########################

library("cmdstanr")
library("kinship2")
library("MBEXCIS")

###################  step 2: viewing the structure of the example data  ##############

head(mixture_data1)
##famid iid fid mid sex type   disease   rs123456
##   1   1   0   0   1    1   2.3184299        0
##   1   2   0   0   2    1  -0.8222242        0
##   1   3   1   2   1    1  -0.1111674        0
##   1   4   1   2   1    1  -1.2965940        0
##   2   1   0   0   1    1   5.0742535        2
##   2   2   0   0   2    1   0.5587576        1

##The first five columns of the data are famid (pedigree Id), iid (individual ID), fid (father ID), mid (mother ID) and sex. The fid and the mid of founders or unrelated individuals are both set to be 0. The numerical codes for sex are 0=unknown, 1=male, 2=female. The sixth column is type, which indicate the data type (type=1 for pedigree data; type=2 for unrelated data). The seventh column "disease" is the trait value. The eighth column "rs123456" is the genotype of target SNP, which is coded as 0, 1 and 2, indicating the number of minor alleles.

head(covariate_mixture1_2)
##famid iid fid mid sex covariate1 covariate2
##   1   1   0   0   1      -0.43          1
##   1   2   0   0   2       0.68          0
##   1   3   1   2   1       1.96          0
##   1   4   1   2   1       1.11          0
##   2   1   0   0   1      -0.35          0
##   2   2   0   0   2       0.67          0

##The first five columns should be consistent with that in the data "mixture_data1". The last two columns are the covariates we want to add. "covariate1" is quantitative and "covariate2" is qualitative.

########################  step 3: explaining the corresponding functions  ################

##Users can use the Bayesian methods to estimate the degree of the skewness of XCI based on the mixtures of pedigree and unrelated data through the function "M_Bayes_XCI".
##When the input variable "prior"="normal", the Bayesian method is MBN; when the input variable "prior"="uniform", the Bayesian method is MBU.
##When the input variable "prior"="customize", users could specify the prior distributions of gamma and other unknown parameters according to their own research background.
##Here let "prior" = "normal"

########################  step 4: running the functions  #######################

##Special note:
##The results may be different for different runs, because of the sampling randomness of the HMC algorithm.
##If the fixed results are wanted, the seed number should be set before running the "M_Bayes_XCI".
##Note that different version of R may lead to different results under the same seed number. The results of the examples are obtained under the R with version 4.1.2.
##Because cmdstanr runs HMC sampling in C language, the stan file in R language needs to be compiled into C language before each run, which may take some extra time.

set.seed(123456)
M_Bayes_XCI(mixture_data=mixture_data1, covariate=covariate_mixture1_2, trait_type="quantitative", trait_missing=NA,
            genotype_missing=NA, covariate_missing=NA, gamma_prior="normal",
            prior_customize=NULL, chains_num=4, parallel_chains=4,
            iter_num=2000, warmup_num=1000, acceptance_rate=0.9, decimal=4)
##results:
##SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
##rs123456         0.3694     0.0078     0.8148 1.0006

######################  step 5: interpreting the results  ######################

##The point estimate of gamma obtained by the MBN method is 0.3694, and the 95% HPDI of gamma is (0.0078, 0.8148).
##The XCI pattern of SNP rs123456 on the trait may be XCI-S towards the minor allele, and where only about 18.47% (0.3694/2) of the cells have the minor allele active overall, and the other 81.53% of the cells have the major allele active overall.
##If the point estimate of gamma is close to 1 and 95% HPDI of gamma contains 1, the XCI pattern of SNP rs123456 on the trait may be XCI-R or XCI-E.
##Rhat is a diagnostic factor assesses convergence of Markov chains. Rhat>1.05 indicates that the Markov chain does not converge, then you should increase "iter_num" and "acceptance_rate", or reconsider prior. In this example, Rhat=1.0006.

######################  Examples 2-8 have the similar steps  ##################

##example 2:
##mixture data for quantitative trait with covariates.
##the prior of gamma is set to uniform distribution and other parameters are set to defaults.
set.seed(123456)
M_Bayes_XCI(mixture_data=mixture_data1,covariate=covariate_mixture1_2,trait_type="quantitative",
            gamma_prior="uniform")
#results:
#SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#rs123456         0.3596      4e-04     0.7724 1.0006

##example 3:
##mixture data for quantitative trait with covariates.
##you want to customize the prior of gamma and other unknown parameters.
##users are required to define the prior of gamma and other parameters by "prior_customize" according to
##their own research background. The format of "prior_customize" must be in stan language, and you can refer
##to https://mc-stan.org/rstan/articles/rstan.html. We give an example as follow:
prior_customize<-"
data {                       // Define the input data type (changes are not recommended)
 int<lower=0> N ;
 int<lower=0> covariate_number ;
 vector[N] y ;
 vector[N] X0 ;
 vector[N] GRM_values ;
 matrix[N,2] X ;
 matrix[N,covariate_number] covariate ;
 }
parameters {                  // Define unknown parameter type (changes are not recommended)
 real<lower=0,upper=2> gamma;
 real beta0 ;
 real beta ;
 vector [covariate_number] beta_covariate;
 real<lower=0> tao ;
 real<lower=0> sigma;
 }
transformed parameters{      // Define the transformed parameters after EVD (changes are not recommended)
 vector[N] mean;
 vector[N] std;
 mean = beta0*X0 + beta*gamma*X[,1] + beta*(2-gamma)*X[,2] + covariate*beta_covariate;
 std = sqrt(tao*GRM_values + sigma^2);
 }
model {                          // Define the prior distribution and likelihood function
 y ~ normal(mean,std);           // likelihood function (changes are not recommended)
 beta0 ~ normal (0, 100);        // prior of intercept (you can customize)
 beta ~ normal (0, 10);          // prior of regression coefficient beta (you can customize)
 beta_covariate ~ normal(0, 10); // prior of regression coefficient of covariates (you can customize)
 gamma ~ normal (1, 2);          // prior of gamma (you can customize)
 tao ~ exponential(1);           // prior of the standard deviation of polygenic effects (you can customize)
 sigma ~ exponential(1);         // prior of the standard deviation of residuals (you can customize)
}"
set.seed(123456)
M_Bayes_XCI(mixture_data=mixture_data1,covariate=covariate_mixture1_2,trait_type="quantitative",
            gamma_prior="customize",prior_customize=prior_customize)
#results:
#SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#rs123456         0.3405      8e-04     0.7925 1.0006

##example 4:
##mixture data for qualitative trait with covariates.
##you want to customize the prior of gamma and other unknown parameters.
##the format of "prior_customize" between quantitative and qualitative traits are little different.
##We give the example for qualitative traits as follow:
prior_customize<-"
data {                           // Define the input data type (changes are not recommended)
 int<lower=0> N ;
 array[N] int y ;
 int<lower=0> covariate_number ;
 matrix[N,2] X;
 matrix[N,N] GRM_female;
 matrix[N,covariate_number] covariate ;
}
transformed data {         // Define the input data after Cholesky (changes are not recommended)
 matrix[N,N] C;
 C = cholesky_decompose(GRM_female);
}
parameters {                     // Define unknown parameter type (changes are not recommended)
 real<lower=0,upper=2> gamma;
 real beta0 ;
 real beta ;
 real<lower=0> tao ;
 vector [covariate_number] beta_covariate;
 vector[N] z;
}
model {                       // Define the prior distribution and likelihood function
 vector[N] p;
 vector[N] d;
 z ~ normal (0,1);
 d = tao *( C * z );
 p = beta0 + beta*gamma*X[,1] + beta*(2-gamma)*X[,2] + covariate*beta_covariate + d;
 y ~ bernoulli_logit(p);         // likelihood function (changes are not recommended)
 beta0 ~ normal(0, 100);         // prior of intercept (you can customize)
 beta ~ normal(0, 10);           // prior of regression coefficient beta (you can customize)
 beta_covariate ~ normal(0, 10); // prior of regression coefficient of covariates (you can customize)
 gamma ~ normal(1, 2);           // prior of gamma (you can customize)
 tao ~ exponential(1);           // prior of the standard deviation of polygenic effects (you can customize)
}"
set.seed(123456)
M_Bayes_XCI(mixture_data=mixture_data3,covariate=covariate_mixture3_4,trait_type="qualitative",
            gamma_prior="customize",prior_customize=prior_customize)
#results:
#SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#rs234567         1.2507     0.7905     1.9196 1.0011

##example 5:
##mixture data for quantitative trait with covariates and missing values.
##the prior of gamma is set to truncated normal distribution and other parameters are set to defaults.
set.seed(123456)
M_Bayes_XCI(mixture_data=mixture_data2,covariate=covariate_mixture1_2,trait_type="quantitative",
            gamma_prior="normal")
#results:
#SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#rs123456         0.4527     0.0079     1.0871 1.0006

##example 6:
##mixture data for qualitative trait with covariates and missing values.
##the prior of gamma is set to truncated normal distribution and other parameters are set to defaults.
set.seed(123456)
M_Bayes_XCI(mixture_data=mixture_data4,covariate=covariate_mixture3_4,trait_type="qualitative",
            gamma_prior="normal")
#results:
#SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#rs234567         1.3081     0.7731     1.9549 1.0005

##example 7:
##mixture data for quantitative trait without covariates.
##the prior of gamma is set to truncated normal distribution and other parameters are set to defaults.
set.seed(123456)
M_Bayes_XCI(mixture_data=mixture_data1,covariate=NULL,trait_type="quantitative",
            gamma_prior="normal")
#results:
#SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#rs123456         0.3941     0.0074     0.8238 1.0007

##example 8:
##mixture data for qualitative trait without covariates.
##the prior of gamma is set to truncated normal distribution and other parameters are set to defaults.
set.seed(123456)
M_Bayes_XCI(mixture_data=mixture_data3,covariate=NULL,trait_type="qualitative",
            gamma_prior="normal")
#results:
#SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#rs234567         1.2463     0.7823     1.8924 1.0006



###########################################################################
##################    examples for P_Bayes_XCI()    #########################
###########################################################################

##example 1:
##There is a data set "pedigree_data1" that includes a quantitative trait and the genotype of target SNP of 30 pedigrees and 260 pedigree-related individuals.
##There is the corresponding covariate data "covariate_pedigree1_2" for above individuals.

############################  step 1: loading the R package  ########################

library("cmdstanr")
library("kinship2")
library("MBEXCIS")

###################  step 2: viewing the structure of the example data  ##############

head(pedigree_data1)
##famid iid fid mid sex    disease   rs123456
##   1   1   0   0   1   2.3184299        0
##   1   2   0   0   2  -0.8222242        0
##   1   3   1   2   1  -0.1111674        0
##   1   4   1   2   1  -1.2965940        0
##   2   1   0   0   1   5.0742535        2
##   2   2   0   0   2   0.5587576        1

##The first five columns of the data are famid (pedigree Id), iid (individual ID), fid (father ID), mid (mother ID) and sex. The fid and the mid of founders or unrelated individuals are both set to be 0. The numerical codes for sex are 0=unknown, 1=male, 2=female. The sixth column "disease" is the trait value. The seventh column "rs123456" is the genotype of target SNP, which is coded as 0, 1 and 2, indicating the number of minor alleles.

head(covariate_pedigree1_2)
## famid iid fid mid sex covariate1 covariate2
##    1   1   0   0   1      -0.49          0
##    1   2   0   0   2       0.32          0
##    1   3   1   2   1       1.46          1
##    1   4   1   2   1       1.54          0
##    2   1   0   0   1      -0.34          1
##    2   2   0   0   2      -1.08          0

##The first five columns should be consistent with that in the data "pedigree_data1". The last two columns are the covariates we want to add. "covariate1" is quantitative and "covariate2" is qualitative.

########################  step 3: explaining the corresponding functions  ################

##Users can use the Bayesian methods to estimate the degree of the skewness of XCI based on only pedigree data through the function "P_Bayes_XCI".
##When the input variable "prior"="normal", the Bayesian method is PBN; when the input variable "prior"="uniform", the Bayesian method is PBU.
##When the input variable "prior"="customize", users could specify the prior distributions of gamma and other unknown parameters according to their own research background.
##Here let "prior" = "normal"

########################  step 4: running the functions  #######################

##Special note:
##The results may be different for different runs, because of the sampling randomness of the HMC algorithm.
##If the fixed results are wanted, the seed number should be set before running the "P_Bayes_XCI".
##Note that different version of R may lead to different results under the same seed number. The results of the examples are obtained under the R with version 4.1.2.
##Because cmdstanr runs HMC sampling in C language, the stan file in R language needs to be compiled into C language before each run, which may take some extra time.

set.seed(123456)
P_Bayes_XCI(pedigree_data=pedigree_data1, covariate=covariate_pedigree1_2, trait_type="quantitative", trait_missing=NA,
            genotype_missing=NA, covariate_missing=NA, gamma_prior="normal",
            prior_customize=NULL, chains_num=4, parallel_chains=4,
            iter_num=2000, warmup_num=1000, acceptance_rate=0.9, decimal=4)
#results:
#SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#rs123456         0.6351          0     1.4811 1.0004

######################  step 5: interpreting the results  ######################

##The point estimate of gamma obtained by the PBN method is 0.6351, and the 95% HPDI of gamma is (0, 1.4811).
##The XCI pattern of SNP rs123456 on the trait may be XCI-R or XCI-E.
##If the point estimate of gamma is significantly different from 1 and 95% HPDI of gamma does not contain 1, the XCI pattern of SNP rs123456 on the trait may be XCI-S.
##Rhat is a diagnostic factor assesses convergence of Markov chains. Rhat>1.05 indicates that the Markov chain does not converge, then you should increase "iter_num" and "acceptance_rate", or reconsider prior. In this example, Rhat=1.0006.

######################  Examples 2-8 have the similar steps  ##################

#' ##example 2:
#' ##pedigree data for quantitative trait with covariates.
#' ##the prior of gamma is set to uniform distribution and other parameters are set to defaults.
#' set.seed(123456)
#' M_Bayes_XCI(pedigree_data=pedigree_data1,covariate=covariate_pedigree1_2,trait_type="quantitative",
#'             gamma_prior="uniform")
#' #results:
#' #SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#' #rs123456         0.5972          0     1.5246 1.0016
#'
#' ##example 3:
#' ##pedigree data for quantitative trait with covariates.
#' ##you want to customize the prior of gamma and other unknown parameters.
#' ##users are required to define the prior of gamma and other parameters by "prior_customize".
#' ##The example of "prior_customize" is same as that in example 3 of function "M_Bayes_XCI".
#' set.seed(123456)
#' P_Bayes_XCI(pedigree_data=pedigree_data1,covariate=covariate_pedigree1_2,trait_type="quantitative",
#'             gamma_prior="customize",prior_customize=prior_customize)
#' #results:
#' #SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#' #rs123456         0.5728     0.0018     1.4738 1.0006
#'
#' ##example 4:
#' ##pedigree data for qualitative trait with covariates.
#' ##you want to customize the prior of gamma and other unknown parameters.
#' ##users are required to define the prior of gamma and other parameters by "prior_customize".
#' ##The example of "prior_customize" is same as that in example 4 of function "M_Bayes_XCI".
#' set.seed(123456)
#' P_Bayes_XCI(pedigree_data=pedigree_data1,covariate=covariate_pedigree1_2,trait_type="qualitative",
#'             gamma_prior="customize",prior_customize=prior_customize)
#' #results:
#' #SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#' #rs234567         1.1429     0.5193     1.9999 1.0001
#'
#' ##example 5:
#' ##pedigree data for quantitative trait with covariates and missing values.
#' ##the prior of gamma is set to truncated normal distribution and other parameters are set to defaults.
#' set.seed(123456)
#' P_Bayes_XCI(pedigree_data=pedigree_data2,covariate=covariate_pedigree1_2,trait_type="quantitative",
#'             gamma_prior="normal")
#' #results:
#' #SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#' #rs123456         0.9163     0.1115     1.9141 1.0001
#'
#' ##example 6:
#' ##pedigree data for qualitative trait with covariates and missing values.
#' ##the prior of gamma is set to truncated normal distribution and other parameters are set to defaults.
#' set.seed(123456)
#' P_Bayes_XCI(pedigree_data=pedigree_data4,covariate=covariate_pedigree3_4,trait_type="qualitative",
#'             gamma_prior="normal")
#' #results:
#' #SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#' #rs234567         1.1992     0.5012     1.9757 1.0004
#'
#' ##example 7:
#' ##pedigree data for quantitative trait without covariates.
#' ##the prior of gamma is set to truncated normal distribution and other parameters are set to defaults.
#' set.seed(123456)
#' P_Bayes_XCI(pedigree_data=pedigree_data1,covariate=NULL,trait_type="quantitative",
#'             gamma_prior="normal")
#' #results:
#' #SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#' #rs123456         0.6802     0.0202     1.4178 1.0011
#'
#' ##example 8:
#' ##pedigree data for qualitative trait without covariates.
#' ##the prior of gamma is set to truncated normal distribution and other parameters are set to defaults.
#' set.seed(123456)
#' P_Bayes_XCI(pedigree_data=pedigree_data3,covariate=NULL,trait_type="qualitative",
#'             gamma_prior="normal")
#' #results:
#' #SNP_Name Point_Estimate HPDI_Lower HPDI_Upper   Rhat
#' #rs234567         1.1201     0.5049     1.9841 1.0021



