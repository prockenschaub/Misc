
# This script adds logistic regression with fractional polynomials 


# Register the model ------------------------------------------------------

set_model_engine("logistic_reg", "classification", eng = "mfp")
set_dependency("logistic_reg", eng = "mfp", pkg = "mfp")


# Define the fit function -------------------------------------------------


mfp_from_dataframe <- function(y, x, df, ...){
  # Parses the standard parameters provided by the ``parsnip`` package to 
  # conform with the formula interface required by ``mfp``.
  #
  # Parameters
  # ----------
  #  y : vector
  #    vector of target values
  #  x : data.frame
  #    predictors (excluding the intercept)
  #  df : scalar integer
  #    degrees of freedom for the fractional polynomial terms
  
  library(mfp)
  
  # Determine which variables should be modelled as fractional polynomials
  vars <- names(x)
  nums <- map_lgl(x, is.numeric)
  terms <- if_else(nums, paste0("fp(", vars, ", df = ", df, ")"), vars)
  str_form <- paste0("y ~ ", paste0(terms, collapse = "+"))
  
  # Call `mfp` with the updated formula and fit the model
  mfp(formula = as.formula(str_form), data = cbind(y = y, x), family = binomial, ...)
}

set_fit(
  model = "logistic_reg",
  eng = "mfp",
  mode = "classification",
  value = list(
    interface = "data.frame",
    protect = c("formula", "data", "family"),
    func = c(fun = "mfp_from_dataframe"),
    defaults = list(df = 2)
  )
)



# Define the predict function ---------------------------------------------

# Probabilities
set_pred(
  model = "logistic_reg",
  eng = "mfp",
  mode = "classification",
  type = "prob",
  value = list(
    pre = NULL,
    post = function(x, object) {
      x <- tibble(v1 = 1 - x, v2 = x)
      colnames(x) <- object$lvl
      x
    },
    func = c(fun = "predict"),
    args =
      list(
        object = quote(object$fit),
        newdata = quote(new_data),
        type = "response"
      )
  )
)


# Classes
set_pred(
  model = "logistic_reg",
  eng = "mfp",
  mode = "classification",
  type = "class",
  value = list(
    pre = NULL,
    post = parsnip:::prob_to_class_2,
    func = c(fun = "predict"),
    args =
      list(
        object = quote(object$fit),
        newdata = quote(new_data),
        type = "response"
      )
  )
)

