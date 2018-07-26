#' Description of the function
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


                          private$df <-self$res %>%
                            select(Beta) %>%
                            map( function(x){
                            x %>% as.matrix() %>%
                              as.data.frame() %>%
                              summarise_all(~sum(.!=0))
                          })
                          },

                          sel_cv = function(s = "lambda.min"){
                            if(is.null(private$cv) || (s != private$s)){
                              private$s <- s
                              self$cv <- private$y %>%
                                gather() %>%
                                group_by(key) %>%
                                nest(.key = Data) %>%
                                mutate(Data = map(Data, ~.$value)) %>%
                                mutate(Model = purrr::map(Data,
                                                          ~ cv.glmnet(private$x, .)))

                            }

                           self$cv <- self$cv %>%
                             mutate(Beta = map(Model, ~coef(.,"lambda.min")[-1])) %>%
                              select(Beta, Data) %>% mutate(Beta = set_names(Beta, colnames(private$y)),
                                                            Data = set_names(Data, colnames(private$y))) %>%

                             summarise_all(~list(bind_rows(.)))

                          }
                        ))
