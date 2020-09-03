
# singcar

<!-- badges: start -->
<!-- badges: end -->

The aim of the R package `singcar` is to provide and encourage usage of
appropriate statistical methods for comparing a case against a control sample.
For instance, they may commonly be done in a neuropsychological context, in
which an individual has incurred a specific brain injury and we wish to test
whether this damage has led to an impairment of some cognitive function and
whether two different functions are dissociable. For many functions there is
normed data available which the patient can be compared against directly.
However, when this is not possible a control sample estimating the population,
against which we wish to compare the patient, must be used. Both frequentist and
Bayesian methods have been developed to do this, first and foremost by John
Crawford and Paul Garthwaite (Crawford et al., 2011; Crawford & Garthwaite,
2002, 2007, 2005; Crawford & Howell, 1998). It is these methods that `singcar`
implements. Power calculators for these tests are also provided. Although the
canonical applications for these tests are in Cognitive Neuropsychology or
Clinical Neuropsychology, they are potentially applicable to any circumstance in
which a measure taken from a single individual is to be compared against data
from a normative sample (i.e. a control group). It should be noted that these
statistical methods could also be applied as a general method of outlier
detection in small samples. 

## Installation

You can install the unstable(!) developmental version of `singcar` by running the following:

```R
install.packages("devtools")
library("devtools")
install_github("jorittmo/singcar")
library("singcar")
```

## Example

The package comes with the dataset `size_weight_illusion`, a neuropsychological
dataset from an investigation of the size-weight illusion in DF, a patient with
visual form agnosia following following bilateral lesions to the lateral
occipital complex (Hassan et al., 2020). It was investigated whether
DF experienced visual size-weight illusion to the same extent as controls (n = 28)
and whether visual and kinesthetic size-weight illusion could be dissociable.
Below follows examples of how to analyse this dataset using the tests provided
in `singcar`.

### Testing for a deficit

If we want to assess whether DF has an impairment compared to controls on
visual size-weight illusion we can test this using a modified two-sample t-test,
called TD (test of deficit: Crawford & Howell, 1998).

``` r
library(singcar)
# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

TD(case = DF_V_SWI, controls = CON_V_SWI, conf_int = TRUE)

##  	Crawford-Howell (1998) t-test

## data:  case = 0.03 and controls (M = 0.16, SD = 0.08, N = 28)
## t = -1.7243, df = 27, p-value = 0.04804
## alternative hypothesis: true difference between case and controls is less than 0
## sample estimates:
## Standardised case score (Z-CC), 95% CI [-2.34, -1.15]       Proportion below case (%), 95% CI [0.95, 12.47] 
##                                             -1.754857                                              4.804003 

```
This can similarily be tested with a Bayesian version of the same test,
yielding approximately (since this test is based on MCMC methods) the same output (Crawford & Garthwaite, 2007).

``` r
# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient 
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

BTD(case = DF_V_SWI, controls = CON_V_SWI)

##  	Bayesian Test of deficit by Crawford and Garthwaite (2007)

## data:  case = 0.03 and controls (M = 0.16, SD = 0.08, N = 28)
## est. z = -1.7388, df = 27, p-value = 0.04803
## alternative hypothesis: true difference between case and controls is less than 0
## sample estimates:
## Std. case score (Z-CC), 95% credible interval [-2.34, -1.15]    Proportion below case (%), 95% credible interval [0.96, 12.44] 
##                                                -1.754857                                                          4.802900

```
If the control sample for a study is not appropriately matched to the case on
variables such as e.g. age or education level it is appropriate to use tests
that account for this by allowing for the inclusion of covariates. Including
theoretically sound covariates is often a good idea. To do this
Crawford et al. (2011) extended their Bayesian verison of the TD. This test
assess the patient on the task of interest by essentially comparing him/her to
the controls with the same score on the covariate.

``` r
# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

# Extracting the coviariate below
DF_age <- size_weight_illusion[size_weight_illusion$PPT == "DF", "YRS"] # Patient
CON_age <- size_weight_illusion[size_weight_illusion$PPT != "DF", "YRS"] # Controls

## BTD_cov(case_task = DF_V_SWI, case_covar = DF_age, control_task = CON_V_SWI, control_covar = CON_age)

##  	Bayesian Test of deficit with Covariates

## data:  case = 0.03 and controls (M = 0.16, SD = 0.08, N = 28)
## est. z = -1.6828, p-value = 0.05386
## alternative hypothesis: true difference between case and controls is less than 0
## sample estimates:
## Std. case difference (Z-CCC), 95% credible interval [-2.31, -1.10]     Proportion below case (%), 95% credible interval [1.05, 13.52] 
##                                                          -1.749556                                                           5.386088 

```
### Testing for a dissociation

If we want to assess whether DF has a dissociation between two functions we can
use a modified paired samples t-test to assess the size of the difference
between the case scores from the two tasks to the distribution of differences
between the tasks in the controls. This can however only be done directly using
the t-distribution if the tasks are measured on the same scale and is called the
unstandardised difference test (UDT: Crawford & Garthwaite, 2005). In the
`size_weight_illusion` dataset it is possible to use this test to whether
patient DF exhibits a dissociation between visual size-weight illusion and
kinesthetic size-weight illusion because the visual and kinaesthetic conditions
are parallel versions of the same task, with different sensory cues. This would
be done as shown below:

``` r
library(singcar)
# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

# Extracting scores from the kinesthetic size-weight illusion from size_weight_illusion 
DF_K_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "K_SWI"] # Patient
CON_K_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "K_SWI"] # Controls

UDT(case_a = DF_V_SWI, case_b = DF_K_SWI, controls_a = CON_V_SWI, controls_b = CON_K_SWI)

##  	Unstandardised Difference Test

## data:  Case score A: 0.03, Case score B: 0.10, Controls A (mean, sd): (0.16, 0.08), Controls B (mean, sd): (0.18, 0.10)
## t = -0.6667, df = 27, p-value = 0.5106
## alternative hypothesis: true difference between tasks is not equal to 0
## sample estimates:
##                       Standardised case score, task A (Z-CC)                       Standardised case score, task B (Z-CC) 
##                                                  -0.13647439                                                  -0.07931545 
## Standardised task discrepancy (Z-DCC), 95% CI [-1.53, -0.59]              Proportion below case (%), 95% CI [6.35, 27.68] 
##                                                  -1.06478887                                                  25.53097678  

```

Most often this is not possible because we wish to estimate abnormality of 
discrepancy on tasks that are not comparable. So otherwise, that
is if the scores must be standardised to be comparable, a statistic that
approximates the t-distribution has been developed and should be used (the
revised standardised difference test RSDT: Crawford & Garthwaite, 2005). 
The visual and kinesthetic size-weight illusion will be used for illustrative
purposes here as well:

``` r
# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

# Extracting scores from the kinesthetic size-weight illusion from size_weight_illusion 
DF_K_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "K_SWI"] # Patient
CON_K_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "K_SWI"] # Controls

RSDT(case_a = DF_V_SWI, case_b = DF_K_SWI, controls_a = CON_V_SWI, controls_b = CON_K_SWI)

##  	Revised Standardised Difference Test

## data:  Case score A: 0.03, Case score B: 0.10, Controls A (mean, sd): (0.16, 0.08), Controls B (mean, sd): (0.18, 0.10)
## approx. abs. t = 1.015, df = 27, p-value = 0.3191
## alternative hypothesis: true difference between tasks is not equal to 0
## sample estimates:
##                         Case score on task A as standard (z) score                         Case score on task B as standard (z) score 
##                                                        -1.7548574                                                         -0.7836956 
##  Std. effect size (Z-DCC) for task diff. between case and controls Proportion of control population with more extreme task difference 
##                                                         -1.0647889                                                         15.9560625 


```
A Bayesian version of this test was also developed [@crawford_comparison_2007],
however, unlike `TD` and `BTD` the `RSDT` and `BSDT` (Bayesian standardised
difference test) differ somewhat and `BSDT` has been shown to keep a better
control of Type I errors if a patient exhibits extreme deficits on both tasks of
interest. Therefore the `BSDT` is recommended above `RSDT`. The usage of the two
R functions is very similar. Since the `BSDT` is based on MCMC methods it can be
quite computationally intensive, depending on the number of iterations you choose.


``` r

# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

# Extracting scores from the kinesthetic size-weight illusion from size_weight_illusion 
DF_K_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "K_SWI"] # Patient
CON_K_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "K_SWI"] # Controls

BSDT(case_a = DF_V_SWI, case_b = DF_K_SWI, controls_a = CON_V_SWI, controls_b = CON_K_SWI, iter = 10^6)

##  	Bayesian Standardised Difference Test

## data:  Case score A: 0.03, Case score B: 0.10, Controls A (mean, sd): (0.16, 0.08), Controls B (mean, sd): (0.18, 0.10)
## est. z = -1.0736, p-value = 0.3067
## alternative hypothesis: true difference between tasks is not equal to 0
## sample estimates:
##                                                              Case score on task A as standard (z) score 
##                                                                                              -1.7548574 
##                                                              Case score on task B as standard (z) score 
##                                                                                              -0.7836956 
## Std. effect size (Z-DCC) for task diff. between case and controls, 95% credible interval [-1.71, -0.45] 
##                                                                                              -1.0647889 
## Proportion of control population with more extreme task difference, 95% credible interval [4.32, 32.47] 
##                                                                                              15.3357743 


```

Just as for `BTD` a version of `BSDT` allowing for
covariates has been developed. This test assess the patient on the discrepancy
between the tasks of interest by essentially comparing him/her to the controls
with the same score on the covariate. 

``` r
# Extracting scores from the visual size-weight illusion from size_weight_illusion 
DF_V_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "V_SWI"] # Patient
CON_V_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "V_SWI"] # Controls

DF_K_SWI <- size_weight_illusion[size_weight_illusion$PPT == "DF", "K_SWI"] # Patient
CON_K_SWI <- size_weight_illusion[size_weight_illusion$PPT != "DF", "K_SWI"] # Controls

# Extracting the coviariate below
DF_age <- size_weight_illusion[size_weight_illusion$PPT == "DF", "YRS"] # Patient
CON_age <- size_weight_illusion[size_weight_illusion$PPT != "DF", "YRS"] # Controls

BSDT_cov(case_tasks = c(DF_V_SWI, DF_K_SWI ), case_covar = DF_age,
         control_tasks = cbind(CON_V_SWI, CON_K_SWI), control_covar = CON_age, iter = 10^5)
         
##  	Bayesian Standardised Difference Test with Covariates

## data:  Case score Y1: 0.03, Case score Y2: 0.10, Controls score Y1: 0.16, Controls score Y2: 0.18
## ave. z = -1.0229, p-value = 0.3303
## alternative hypothesis: true difference between tasks is not equal to 0
## sample estimates:
##                                                               Case score on task X as standard (z) score 
##                                                                                                -1.754857 
##                                                               Case score on task Y as standard (z) score 
##                                                                                                -0.783696 
## Std. effect size (Z-DCCC) for task diff. between case and controls, 95% credible interval [-1.67, -0.40] 
##                                                                                                -1.064152 
##  Proportion of control population with more extreme task difference, 95% credible interval [4.77, 34.54] 
##                                                                                                16.510000       

```
All of the functions above can also take summary (mean, sd, control sample size) data as input.

### Power calculators

A further capacity of `singcar` is that it can be used to calculate power for
for these single case-control comparisons. Calculations for all Bayesian tests
and `RSDT` are simulation based and (especially the tests with covariates) can
be computationally intense. Calculators for `TD` and `UDT` (unstandardised
difference test) are exact (their power functions have been derived
analytically) and can both be used to find a specific sample size given a
desired power. For the other calculators all parameters must be given. Means and
standard deviations for the control population are at default set to 0 and 1
meaning that the case value will be interpreted as differences from the mean in
standard deviations, these parameter values can be changed as you like. Examples
are given below:

``` r
TD_power(case = -2, power = 0.8, mean = 0, sd = 1, alternative = "two.sided")

## [1] "Power(0.44280) will not increase more than 0.5% for any additional participant over n = 16"
##    n     power
## 1 16 0.4428042

TD_power(case = 70, sample_size = 10, mean = 100, sd = 15, alternative = "less", alpha = 0.1)

## [1] 0.7039033

RSDT_power(case_a = 70, case_b = 20, mean_a = 100, mean_b = 25, sd_a = 15, sd_b = 10, sample_size = 10) 

## [1] 0.5689


# Takes long time to compute
BTD_cov_power(case = -2, case_cov = 0, control_task = c(0, 1),
              control_covar = c(0, 1), cor_mat = diag(2), sample_size = 10)

# [1] 0.511

```

# References

Crawford, J., & Garthwaite, P. (2002). Investigation of the single case in neuropsychology: Confidence limits on the abnormality of test scores and test score differences. Neuropsychologia, 40(8), 1196–1208. https://doi.org/10.1016/S0028-3932(01)00224-X

Crawford, J., & Garthwaite, P. (2007). Comparison of a single case to a control or normative sample in neuropsychology: Development of a Bayesian approach. Cognitive Neuropsychology, 24(4), 343–372. https://doi.org/10.1080/02643290701290146

Crawford, J., & Garthwaite, P. (2005). Testing for Suspected Impairments and Dissociations in Single-Case Studies in Neuropsychology: Evaluation of Alternatives Using Monte Carlo Simulations and Revised Tests for Dissociations. Neuropsychology, 19(3), 318–331. https://doi.org/10.1037/0894-4105.19.3.318

Crawford, J., Garthwaite, P., & Ryan, K. (2011). Comparing a single case to a control sample: Testing for neuropsychological deficits and dissociations in the presence of covariates. Cortex, 47(10), 1166–1178. https://doi.org/10.1016/j.cortex.2011.02.017

Crawford, J., & Howell, D. (1998). Comparing an Individual’s Test Score Against Norms Derived from Small Samples. The Clinical Neuropsychologist, 12(4), 482–486. https://doi.org/10.1076/clin.12.4.482.7241

Hassan, E. K., Sedda, A., Buckingham, G., & McIntosh, R. D. (2020). The size-weight illusion in visual form agnosic patient DF. Neurocase, 1–8. https://doi.org/10.1080/13554794.2020.1800748
