#' Load all packages needed, installing any missing packages from CRAN.
#'
#' @param auto_install Install any packages that could not be loaded.
#' @param update Update any packages that can be updated.
#' @param java_mem Amount of RAM to allocate to rJava; must happen before
#'   library is loaded.
load_all_packages = function(auto_install = F, update = F, java_mem = "4g", verbose = F) {
  if (verbose) {
    # Output R version so we know which package versions we're using.
    cat(R.version.string, "\n")
  }

  # Allocate 4GB to Java for bartMachine; needs to happen before we load rJava library.
  options(java.parameters = paste0("-Xmx", java_mem))

  libs = c(
    # Warning: This will crash R if a JDK is not installed.
    "bartMachine",
    "caret",
    "doMC",
    #"doParallel",
    #"doSNOW",
    #"e1071",
    "earth",
    #"gam",
    #"gbm",
    #"gee",
    "ggplot2",
    "glmnet",
    #"KernelKnn",
    "kernlab",
    "mgcv",
    "nnet",
    #"pROC",
    "randomForest",
    "ranger",
    "RhpcBLASctl",
    #"ROCR",
    #"rpart",
    #"sandwich",
    #"speedglm",
    "xgboost"#,
    #"xtable"
  )

  # Code is not yet run. We run afterward, possibly with messages suppressed.
  expression = quote({

    # Install devtools if we don't already have it.
    if (!require("devtools") && auto_install) {
      if (verbose) cat("No devtools detected - installing.\n")
      install.packages("devtools")
      library(devtools)
    }


    # NOTE: may want to install the latest xgboost from github.
    # Can run this manually:
    if (!require("xgboost") && auto_install) {
      if (verbose) cat("No xgboost detected - installing.\n")
      install.packages("xgboost",
                       repos=c("http://dmlc.ml/drat/", getOption("repos")),
                       type="source")
    }

    # Install any packages from github that are needed.
    if (auto_install) {
      devtools::install_github(c("ecpolley/SuperLearner",
                                 "ck37/ck37r"))
    }
    invisible(sapply(c("SuperLearner", "ck37r"), require, character.only = T))

    ck37r::load_packages(libs, auto_install, update, verbose = verbose)
  })

  # Now run the stored code either directly or with messages suppressed.
  if (verbose) {
    # Allow messages to be output.
    eval(expression)
  } else {
    # Supress messages.
    suppressMessages(eval(expression))
  }
}
