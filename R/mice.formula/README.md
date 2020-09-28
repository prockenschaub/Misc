# Comparison of `mice`'s formula and predictorMatrix interface for interactions

A recent [Master's project](https://uknowledge.uky.edu/cgi/viewcontent.cgi?article=1234&context=cph_etds) reported different results when specifying an imputation model with interaction using `mice`'s `predictorMatrix` and `formulas` arguments, claiming that using `predictorMatrix` results in more precise estimates. This is surprising because in theory both options were designed to result in equivalent models. The code in this folder revisits the code published as part of the Master's thesis and investigates the reasons for a difference in results.

## Findings 

Results were close to those reported in the report could be reproduced using the original code (unfortunately no seeds were set, so results couldn't be reproduced exactly). Differences in estimates between `predictMatrix` and `formulas` could be traced back to an incorrect definition of the `predictorMatrix`, which included an interaction between `x` and `y` (`zyint`) in the imputation of `y`.


```r
myPredMat
#>        x y.miss z.miss xzint yzint xyint
#> x      0      1      1     1     1     1
#> y.miss 1      0      1     1     0     1
#> z.miss 1      1      0     0     0     1
#> xzint  1      1      1     0     1     1
#> yzint  1      1      1     1     0     1
#> xyint  1      1      1     1     1     0
```

After correcting the `predictMatrix`, both interfaces provided identical results.



## File structure

- formula.R: original code using the `formulas` interface (pp. 13-17)
- predict_mat.R: original code using the `predictMatrix` interface (pp. 26-31)
- predict_mat2.R: corrected code using the `predictMatrix` interface 
- graph.R: original code to create plots (pp. 41-50)

Folders store the outcome of the simulations.