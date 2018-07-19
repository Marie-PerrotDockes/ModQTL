#' Description of the function
#'
#' @param  response a vector response variable
#' @param  regressors a quantitative matrix of regressor
#' @param  group a vector with two levels. (The group of the ANCOVA)
#' @param  a the parameters that indicate how much the coefficients will be fused
#' @param  lambda if the user wants to use it owns values of lambdas
#' @return The coefficients of the fused lasso ANCOVA for the different value of lambda
#' @examples
#' B <- c(1, -1, 1.5, 1.5, rep(0, 6), 2, 0, 2, 0)
#'group <- c(rep('M1', 10), rep('M2', 10))
#'regressors <- matrix(rnorm(6*20), ncol = 6)
#'X  <- model.matrix(~group + group:regressors - 1)
#'y <- X%*%B + rnorm(20)
#'y <- scale(y)
#'mod <- fl2(y, regressors, group)
#'colors <- c(rep("grey",2), rep('green',2),rep('black', 6), rep(c("orange","blue"), 2), 'darkgreen', rep('yellow',3), rep('purple',2))
#'matplot(mod$lambda ,t(mod$beta),type='l',col=colors)
#' @export
QTLmod_lasso <- R6Class("QTLmod_lasso",
                        inherit = QTLmod,
                        public = list(
                          initialize = function(x, y){
                            super$initialize(x, y)
                          },
                          estime = function(lambda = NULL, standardize = FALSE, intercept = TRUE){
                          self$mod <-  private$y %>%
                              gather() %>%
                              group_by(key) %>%
                              nest(.key = Data) %>%
                              mutate(Data = map(Data, ~.$value)) %>%
                              mutate(Model = purrr::map(Data,
                                                        ~ glmnet(private$x, .,
                                                                 lambda = lambda,
                                                                 standardize = standardize,
                                                                 intercept = intercept)))
                            self$res <- self$mod  %>%
                              mutate( Beta = map(Model, ~.$beta), Intercept = map(Model,~.$a0)) %>%
                              select(-Model)
                                                    }

                        ))
