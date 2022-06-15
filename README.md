
<!-- README.md is generated from README.Rmd. Please edit that file -->

 Please be aware that the RNAseqStat2 https://github.com/xiayh17/RNAseqStat2 is an new version of RNAseqStat.
 RNAseqStat won't be updated any more.

# RNAseqStat

<!-- badges: start -->

[![codecov](https://codecov.io/gh/xiayh17/RNAseqStat/branch/master/graph/badge.svg?token=6rlirZDUVo)](https://codecov.io/gh/xiayh17/RNAseqStat)

<!-- badges: end -->

The goal of RNAseqStat is a workflow for DEG analysis.

![RNAseqStat](https://cdn.jsdelivr.net/gh/xiayh17/Figs@main/uPic/RNAseqStat.svg)

## Installation

You can install the released version of RNAseqStat from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("RNAseqStat")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xiayh17/RNAseqStat")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(RNAseqStat)
#> 
## basic example code
```

## run all in one time

`runAll` will run full workflow in default parameters.

This include 5 steps:

1.  Check your data
2.  DEG analysis
3.  GO
4.  KEGG
5.  GSEA

``` r
runAll(count_data = counts_input, group_list = group_list, 
       test_group = "T", control_group = "C",
       OrgDb = 'org.Hs.eg.db', dir = results_dir)
```

## Every step can be run separately.

Your can check more in xiayh17.top/rnaseqstat

![paste-912DEDCA](https://cdn.jsdelivr.net/gh/xiayh17/Figs@main/uPic/paste-912DEDCA.png)

![paste-AED756ED](https://cdn.jsdelivr.net/gh/xiayh17/Figs@main/uPic/paste-AED756ED.png)
