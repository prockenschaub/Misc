
# This script adds logistic regression with highly adaptive lasso


# Register the model ------------------------------------------------------

set_model_engine("logistic_reg", "classification", eng = "hal")
set_dependency("logistic_reg", eng = "hal", pkg = "hal9001")

set_model_arg(
  model = "logistic_reg",
  eng = "hal",
  parsnip = "penalty",
  original = "lambda",
  func = list(pkg = "dials", fun = "penalty"),
  has_submodel = FALSE
)


# Define the fit function -------------------------------------------------

port_hal <- function(x, y, ...){
  # Parses the standard parameters provided by the ``parsnip`` package to 
  # conform with the parameter names used in ``hal9001``.
  #
  # Parameters
  # ----------
  #  y : vector
  #    vector of target values
  #  x : data.frame
  #    predictors (excluding the intercept)
  
  library(hal9001)
  fit_hal(X = x, Y = y, ...)
}

set_fit(
  model = "logistic_reg",
  eng = "hal",
  mode = "classification",
  value = list(
    interface = "matrix",
    protect = c("x", "y", "X", "Y"),
    func = c(fun = "port_hal"),
    defaults = list(family = "binomial",
                    X_unpenalized = NULL,
                    cv_select = FALSE,
                    yolo = FALSE)
  )
)


# Define the predict function ---------------------------------------------

# Probabilities
set_pred(
  model = "logistic_reg",
  eng = "hal",
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
        new_data = quote(new_data),
        type = "response"
      )
  )
)


# Classes
set_pred(
  model = "logistic_reg",
  eng = "hal",
  mode = "classification",
  type = "class",
  value = list(
    pre = NULL,
    post = parsnip:::prob_to_class_2,
    func = c(fun = "predict"),
    args =
      list(
        object = quote(object$fit),
        new_data = quote(new_data),
        type = "response"
      )
  )
)
