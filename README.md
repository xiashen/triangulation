
## Installation

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
- The methods are **conditionally independent**, i.e., given a single test result of one method, we know nothing about the result of another method.
