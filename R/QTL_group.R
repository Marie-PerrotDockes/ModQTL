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

                            self$reord_beta()



                            private$df <- self$res %>%
                              select(Beta) %>%
                              map( function(x){
                                x %>%
                                  extract2(1) %>%
                                  as.data.frame() %>%
                                  group_by(self$group[1:ncol(private$x)]) %>%
                                  summarise_all(sum) %>%
                                  summarise_all(~sum(.!=0))
                              }) %>%
                              map(function(x){
                                x %>%
                                  select(-1) %>%
                                  as.numeric()
                              })
                          },
                          sel_cv = function( s = "lambda.min"){
                            self$cv <-  private$tb %>%
                              mutate(Model = map(Data,
                                                 ~cv_grp_lasso(private$x, .,
                                                               self$group)
                              )) %>%
                              mutate(Beta = map(Model, ~coef(.,s)[-1,])) %>%
                              mutate(Beta =map(Beta, function(beta){
                                beta %>% as.data.frame() %>%
                                  rownames_to_column() %>%
                                  separate(rowname, c("Trait", "Marker"), "_", extra ="merge") %>%
                                  spread(Trait ,  ".") %>%
                                  arrange(match( Marker,colnames(private$x)))

                              }))

                          }
                        ))

#' Description of the function
#'
#' @export
QTLmod_group_univ <- R6Class("QTLmod_group_univ", inherit = QTLmod_group,
                             public = list(
                               initialize = function(x,y){
                                 super$initialize(x, y)
                                 self$group <- get_group_marker(private$x, r = 1)

                                 private$tb <- private$y %>%
                                   univ_y()


                               },
                               sel_cv = function(s= "lambda.min"){
                                 super$sel_cv()
                                 self$cv <- self$cv %>%  mutate(Beta = map2(Beta, Trait, function(x,y){
                                   x %>% mutate(key = y ) %>% set_names(c("Marker","value","key"))
                                 }) ) %>%
                                   select(Data, Beta) %>%
                                   mutate(Beta = set_names(Beta, colnames(private$y)),
                                                Data = set_names(Data, colnames(private$y))) %>%
                                   summarise_all(~list(bind_rows(.)))
                               },
                               reord_beta= function(){
                                 self$res <- self$res %>%
                                   mutate(Beta = map (Beta,function(x){
                                     ord_beta(beta =x,name=colnames(private$x))
                                   }))
                               },
                               reord_beta_cv = function(){
                                 self$cv <- self$cv %>% select(Data, Beta) %>%
                                   mutate(Beta = set_names(Beta, colnames(private$y)),
                                          Data = set_names(Data, colnames(private$y))) %>%
                                   summarise_all(~list(bind_rows(.)))

                               }
                             ))
#' Description of the function
#' @export
QTLmod_group_multi <- R6Class("QTLmod_group_multi", inherit = QTLmod_group,
                              private= list(r = NULL),
                              public = list(
                                initialize = function(x, y){
                                  super$initialize(x, y)
                                  private$r <- ncol(private$y)
                                  private$tb <- private$y %>%
                                    list() %>%
                                    tibble(Data = .)
                                },

                                reord_beta = function(){
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
                                },

                                sel_cv = function(s ="lambda.min"){
                                  super$sel_cv(s=s)
                                  self$cv <- self$cv %>%
                                    mutate(Beta = map(Beta,~gather(.,key, value, -Marker)))

                                }
                              ))



#' Description of the function
#'
#' @export
QTLmod_group_multi_both <- R6Class("QTLmod_group_multi_both", inherit = QTLmod_group_multi,
                                   public = list(
                                     initialize = function(x, y){
                                       super$initialize(x, y)
                                       self$group <- get_group_both(private$x, r = private$r)
                                     }
                                   ))

#' Description of the function
#'
#' @export
QTLmod_group_multi_marker <- R6Class("QTLmod_group_multi_marker", inherit = QTLmod_group_multi,
                                     public = list(
                                       initialize = function(x, y){
                                         super$initialize(x, y)
                                         self$group <- get_group_marker(private$x, r = private$r)
                                       }
                                     ))
