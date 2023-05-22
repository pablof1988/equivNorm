# equivNorm
This R package implements methods to prove normality (or even another distributions) using an equivalence approach. Simulation functions can be used to measure the effectiveness (in terms of Type One Error Probability and Power) of pretest to prove normality when testing means difference, in addtition the robustness of t-Student and Wilcoxon text can be obtained

The main purpose of this package is show the functions who computed the results of the paper __Goodness and Lack of fit to pretest normality when means comparison tests__


## Installation instructions
equivNorm package has to be installed with a working R version (>=4.2.0). Installation could take a few minutes on a regular desktop or laptop. Package can be installed from `devtools` package, then it needs to be loaded using `library(equivNorm)`


To install from Github

```{r}
devtools::install_github("pablof1988/equivNorm")
```

## Basic ideas and concepts
equivNorm package provides the following functions:

- **lackFitTest:** Prove if a sample comes from a target distribution through  an equivalence approach. 
- **normequiv:** Prove if a sample comes from a normal distribution through an equivalence approach.
- **rnonorm:** Generates random numbers from non normal distributions using the Fleishman's coefficients.
- **epsilon:** Compute a non arbitrary epsilon value of irrelevance for an equivalence normality test of lack of fit
- **rejectH0:** Simulation function to estimate the Type I error probability and power of hypothesis tests.
- **bothTests & tiepT:** These are complementary simulation function whose main purpose is to help to `rejectH0` function.

## Contribution Guidelines
Contributions are welcome, if you wish to contribute or give ideas to improve the package, please you can contact with maintainer (Pablo Flores) to the addres `p_flores@espoch.edu.ec`, and we can discuss your suggestion.
