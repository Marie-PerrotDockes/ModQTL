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
#'@import R6 Matrix gglasso tidyverse glmnet stabs magrittr gganimate viridis
#' @export
QTLmod <- R6Class("QTLmod",
                  private = list(
                    x = NULL,
                    y = NULL,
                    intercept= NULL,
                    r= NULL,
                    Data = NULL,
                     tresh = NULL),
                  public = list(
                    res = NULL,
                    pe = NULL,
                    mod = NULL,
                    initialize = function(x, y){
                      private$x <- as.matrix(x)
                      private$y <- as.data.frame(y)
                      if(nrow(private$x) != nrow(private$y)) stop("X and Y must have the same number of rows!")
                      private$r <- ncol(private$y)
                    },
                    plot_error = function(print = TRUE){
                      if(is.null(self$res$MSE))
                        self$predict()
                      MSE <- self$res$MSE %>%
                        map_dfc(function(x) x %>% t() %>%  as.data.frame()) %>%
                        set_names(colnames(private$y)) %>%
                        rowid_to_column() %>%
                        gather(key, value,-rowid)

                     pe <- ggplot(data = MSE, aes (x = rowid, color = key, y = value )) +
                       geom_point() +
                       geom_line()+
                       theme_bw()

                     pe

                    },
                    plot_coef = function(tresh, sel = colnames(private$y)){
                      if(is.null(self$mod)) self$estime()

                      if(is.null(private$Data) || private$tresh != tresh){
                      private$Data <- self$res$Beta %>% map(~.[,tresh]) %>%
                        set_names(colnames(private$y)) %>%
                        bind_cols(id = colnames(private$x))  %>%
                        set_names(c(colnames(private$y),"id"))
                      private$tresh <- tresh
                      }

                      private$Data <- private$Data %>%
                        filter_at(vars(sel), any_vars(. != 0)) %>%
                        gather(key, value, -id)

                      pc <- ggplot(data = private$Data, aes( x = id, y = key, fill = value)) +
                        scale_fill_viridis() +
                        geom_tile() +
                        theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ylab("") + xlab("")
                      pc
                    },

                    ROC = function(b){
                      self$res$Beta %>%
                        map(function(x){x  %>%
                            as.matrix() %>%
                            as.data.frame() %>%
                            gather() %>%
                            group_by(key) %>%
                            summarise_all(  funs(TPR(., b = b), FPR(., b = b)))
                          })

                    },
                    plot_anime = function(){

                    },
                    # estime = function(){
                    #   self$res <- self$mod  %>%
                    #     mutate( Beta = map(Model, ~.$beta)) %>%
                    #     select(-Model)
                    # },
                    predict = function(){
                         if(is.null(self$mod))
                     self$estime()
                      self $res <- self$res %>%
                        mutate(Yhat = map2(Beta, Intercept, ~ cbind(1,private$x) %*%
                                             rbind(.y,as.matrix(.x)))) %>%
                        mutate(Ehat = map2(Data, Yhat, ~.x - .y)) %>%
                        mutate(Ehat = map(Ehat, ~as.data.frame(.))) %>%
                        mutate(MSE = map(Ehat, ~summarise_all(.,~sum(.^2)))) %>%
                        mutate(DF = map(Beta, function(x){
                          x %>% as.matrix() %>%
                            as.data.frame() %>%
                            summarise_all(~sum(.!=0))
                        }))
                    }
                  )
)

