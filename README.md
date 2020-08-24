
## Installation

You need **R** (version >= 3.1) on your computer. If you have not installed **R**, you can easily download and install it from https://cran.r-project.org/.

We suggest installing the **triangulation** package using the `install_github` function in the **devtools** package:

```{r}
require(devtools)
install_github('https://github.com/xiashen/triangulation')
```

then load the package via:

```{r}
require(triangulation)
```

## Maximum likelihood estimation of operating characteristics

Without a gold standard, operating characteristics of different testing methods can be estimated by maximum likelihood. As long as:

- **At least 3** distinct methods are applied on the same set of tests;
- The methods are **conditionally independent**, i.e., given a single test result of one method, we know almost nothing about the result of another method.

Here we use a 3-method example to demonstrate how the estimation can be done using the **triangulation** package. The three methods were applied on 4,800 independent binary test problems. Let us denote positive and negative test results as 1 and 0, respectively. Based on the 0-1 pattern of the results by the three methods, the 4,800 binary test problems ended up the data as follows:

```{r}
ny <- c(1051, 179, 1028, 154, 1040, 159, 981, 208)
names(ny) <- c('000', '001', '010', '011', '100', '101', '110', '111')
ny

##  000  001  010  011  100  101  110  111
## 1051  179 1028  154 1040  159  981  208
```

where for 1,051 test problems, all the three methods gave negative answers, and so on. Based on such data, we can quickly estimate the **sensitivity** and **specificity** of each method, and the **prevalence** of true positives among the 4,800 tests:

```{r}
triangulate(counts = ny, ntest = 3)

## Estimation done.

## $sensitivity.est
## [1] 0.46523251 0.46148047 0.11314916

## $specificity.est
## [1] 0.39683161 0.39360390 0.74392541

## $prevalence.est
## [1] 0.76868254
```

We can see that the 1st and 2nd methods have similar performance, while the 3rd method has limited sensitivity but relatively high specificity. The **triangulate()** function also provides the possibility to obtain standard errors of these estimates via `B` bootstrap samples. When `B = 10`:

```{r}
triangulate(counts = ny, ntest = 3, B = 10)

## Estimation done.
## Boostrap standard errors:
## Progress: 100%
## $sensitivity.est
## [1] 0.46523251 0.46148047 0.11314916
## 
## $sensitivity.se
## [1] 0.020265770 0.028815577 0.041079261
## 
## $specificity.est
## [1] 0.39683161 0.39360390 0.74392541
## 
## $specificity.se
## [1] 0.12409867 0.10303815 0.13525799
## 
## $prevalence.est
## [1] 0.76868254
## 
## $prevalence.se
## [1] 0.11471999
```

## Combined FDR from different methods

In practice, although it is essential to assess the operating characteristics of different methods in the absense of a gold standard, we normally want to **meta-analyze** the results from different methods to have more confident conclusions. The **combined.fdr()** function achieves this by reporting a combined FDR assessment for each category of test result.

Based on the estimated **sensitivity**, **specificity**, and **prevalence** above, we have:

```{r}
boot <- triangulate(counts = ny, ntest = 3, B = 10)

## Estimation done.
## Boostrap standard errors:
## Progress: 100%

combined.fdr(boot$sensitivity.est, boot$specificity.est, boot$prevalence.est)

##        000        001        010        011        100        101        110        111
## 0.12058012 0.26802765 0.19822130 0.39600078 0.19299140 0.39370497 0.30008041 0.53867587
```

## Citation
If you use the HDL software, please cite

Yang Z, Xu W, Zhai R, Li T, Ning Z, Pawitan Y, Shen X (2020). Triangulation of analysis strategies links complex traits to specific tissues and cell types.  _Submitted_.

## For Help
For direct R documentation of the **triangulate()** and  **combined.fdr()** functions, you can use a question mark in R, e.g.,

```{r}
?triangulate
```

If you want further discussion or still have questions, please feel free to email the maintainer of **triangulation** via xia.shen@ed.ac.uk.








