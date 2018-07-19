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
QTLmod_group <- R6Class("QTLmod_group",
                        inherit = QTLmod,
                        private = list(
                          tb = NULL
                        ),
                        public = list(
                          group = NULL,
                          estime = function(){
                            self$mod <-  private$tb %>%
                              mutate(Model = map(Data,
                                                 ~grp_lasso(private$x, .,
                                                            self$group)
                              ))
                            self$res <- self$mod  %>%
                                  mutate( Beta = map(Model, ~.$beta),
                                          Intercept = map(Model, ~.$b0) ) %>%
                                  select(-Model)
                          }

                        ))

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
QTLmod_group_univ <- R6Class("QTLmod_group_univ", inherit = QTLmod_group,
                             public = list(
                               estime = function(){
                               self$group <- get_group_marker(X, r = 1)
                                 private$tb <- private$y %>%
                                   univ_y()
                                 super$estime()

                                 self$res
                               }
                             ))
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
QTLmod_group_multi <- R6Class("QTLmod_group_multi", inherit = QTLmod_group,
                              private= list(r = NULL),
                              public = list(
                                initialize = function(x, y){
                                  super$initialize(x, y)
                                  private$r <- ncol(private$y)
                                },
                                estime = function(type = "both"){
                                  if( type == "both")
                                   self$group <- get_group_both(private$x, r = private$r)
                                  if(type == "marker")
                                    self$group <- get_group_marker(private$x, r = private$r)
                                  private$tb <- private$y %>%
                                    list() %>%
                                    tibble(Data = .)
                                  super$estime()
                                  # self$res <- self$res %>%
                                  #   mutate(Beta = map2(Beta, Intercept,
                                  #                      ~trans_beta(.x, .y,
                                  #                                  name = colnames(private$x)))) %>%
                                  #   select(Beta) %>%
                                  #   extract2(1)  %>%
                                  #   bind_cols(univ_y(private$y), . )

                                  #
                                  inter <- self$res %>%
                                    select(Intercept) %>%
                                    rep(.,private$r) %>%
                                    set_names(1:private$r) %>%
                                    map(~extract2(.,1)) %>%
                                    tibble() %>% set_names("Intercept")

                                  self$res <- self$res %>% select(Beta) %>%
                                    map(~trans_beta(., name = colnames(private$x))) %>%
                                    extract2(1) %>%
                                    bind_cols(univ_y(private$y), ., inter)

                                }
                              ))

