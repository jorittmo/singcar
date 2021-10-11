---
title: '`singcar`: Comparing single cases to small samples in `R`'
tags:
  - R
  - Case control comparison
  - Neuropsychology
  - Brain lesion studies
  - Small samples
authors:
  - name: Jonathan Ö. Rittmo
    orcid: 0000-0001-5075-0166
    affiliation: 1
  - name: Robert D. McIntosh
    orcid: 0000-0002-7615-6699
    affiliation: 1
affiliations:
 - name: Human Cognitive Neuroscience, Psychology, University of Edinburgh, UK
   index: 1
date: 11 October 2021
bibliography: paper.bib

---

# Summary

Case-control comparisons are a class of statistical tests allowing researchers
to compare single cases to populations estimated from a sample. Such tests have
wide potential utility, but historically have been applied mostly in the fields
of cognitive and clinical neuropsychology, to infer whether individuals have
suffered significant cognitive changes as the consequence of a brain lesion. One
may wish to estimate whether that individual has abnormally low performance on
some cognitive ability, or if one cognitive ability is abnormally discrepant
with respect to another cognitive ability. John Crawford, Paul Garthwaite and
colleagues have developed several related methods to statistically test for
abnormality on a single variate and abnormality of the difference between two
variates when a single case is compared to a small sample, while controlling the
Type I error rate [e.g.,@crawfordComparingIndividualTest1998;
@crawfordComparingSingleCase2011; @crawfordComparisonSingleCase2007;
@crawfordInvestigationSingleCase2002; @crawfordTestingSuspectedImpairments2005].
This paper presents the `R` package `singcar` in which they are implemented. Due to recent discussion
on the fundamental power limits of these tests [@mcintoshPowerCalculationsSinglecase2020] the package also includes 
associated power calculators. 


# Statement of need

There are many reasons why researchers and clinicians might want to look at
single cases instead of at the average of some group. In certain fields, such as
neuropsychology, this need arises because the pattern of naturally-occurring
brain damage will be unique in each individual case. From a theoretical
perspective, this means that a single patient might be the only available source
of data for a given phenomenon. 
From a practical, clinical perspective,
diagnosis and description of the pattern of cognitive impairment is done at the
individual level. Individual brain-damaged patients are thus often compared to
the healthy population to assess changes in cognitive functioning. If we want to
assess the patient score on some variate Y, for which we do not know the
population parameters, these must be estimated from a sample. Thus, the
single-case of interest is compared to a control sample. There are many other
areas where the application of such methods could also be useful, for example
studies of uncommon human expertise.

As it represents the canonical field for the application of these methods, the
nomenclature of neuropsychology is adopted. An abnormally low score on
a single variate is referred to as a *deficit*, an important concept
for clinical and basic neuropsychology alike. For the latter area another
concept is also considered to be of cardinal importance: the ability to test for
an abnormally large discrepancy between two variates. This is referred to as a
*dissociation*, which is taken to provide evidence for some degree of
functional independence between two cognitive abilities. By charting
dissociations, a cognitive architecture of the mind can be theorized
[@shalliceNeuropsychologyMentalStructure1988].

During the last 20 years, a class of related methods have been developed for
case-control comparisons, allowing researchers to estimate abnormality and test
for deficits and dissociations in the single case, while controlling the Type I
error rate. These tests have been developed mainly by John Crawford and Paul
Garthwaite [e.g.,@crawfordComparingIndividualTest1998;
@crawfordComparingSingleCase2011; @crawfordComparisonSingleCase2007;
@crawfordInvestigationSingleCase2002; @crawfordTestingSuspectedImpairments2005].
John Crawford has provided free software packages to perform these tests, making
them available at https://homepages.abdn.ac.uk/j.crawford/pages/dept/psychom.htm. However,
these are available only as standalone compiled computer programs for Windows
operating systems. Many of these programs require manual input of summary
statistics, and output a static text file and for thorough documentation one must
consult the original publications. 

Our aim is to encourage and simplify usage of these methods by implementing them
in the package `singcar` for the `R` environment
[@rcoreteamLanguageEnvironmentStatistical2020], bringing them together in a
fully documented package with open source code that works across platforms.
Further advantages of `singcar` include an API that has more modifiable test
parameters. It is also possible to automate these tests if multiple analyses need to
be run for the purposes of data analysis or simulation studies. 
The development of Crawford and Garthwaite's methods has been focused around limiting Type I
errors, but to emphasise the importance of considering Type II errors we also provide power
calculators for each test function. Our hope in doing so is to increase awareness of power
for this methodology as well as to aid in the planning and design of experiments 
[@mcintoshPowerCalculationsSinglecase2020].

Note that the `R` package `singlecase` [@matthieusinglecase2008] contains some overlapping
functionality with `singcar`, but it has not been maintained since 2008 and lacks
core functionality such as tests allowing for the inclusion of covariates, and
power calculators. A recent study by @mitchell2020peripheral investigating 
peripheral reaching in Alzheimer’s disease and mild cognitive impairment exemplifies
the uses of these novel functionalities in `singcar`. 

# Functionality

`singcar` contains seven functions to estimate a
case's abnormality compared to a normal population estimated from a small
sample, three of them with regards to a single variate and four with regards to
the discrepancy between two variates. Both frequentist and Bayesian methods are
provided, all developed originally by Crawford and colleagues
[@crawfordComparingIndividualTest1998;
@crawfordComparingSingleCase2011; @crawfordComparisonSingleCase2007;
@crawfordInvestigationSingleCase2002; @crawfordTestingSuspectedImpairments2005].
Of special note for psychological research are the methods allowing the
inclusion of covariates [@crawfordComparingSingleCase2011] using Bayesian
regression techniques. These methods make matching the control sample to the
case less cumbersome. For rationale as well as mathematical and contextual
background of the methods consult the package vignette.

# References
