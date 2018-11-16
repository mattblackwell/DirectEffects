
#' Coefficient Estimates Across Bootstrapped Samples
#'
#'
#' @param seqg A fitted sequential_g estimate, an output of sequential_g.
#' @param boots The number of bootstrap replicates. Defaults to 1000.
#' 
#' @return An object of type \env{seqgboots} which is a matrix with \env{boots} 
#' rows and columns for each coefficient in the \env{seqg} model.
#' 
#' @examples 
#' 
#' data(ploughs)
#' form <- women_politics ~ plow +
#'  agricultural_suitability + tropical_climate + large_animals + rugged | 
#'  years_civil_conflict + years_interstate_conflict  + oil_pc +
#'  european_descent + communist_dummy + polity2_2000 | 
#'  centered_ln_inc + centered_ln_incsq
#' s1 <- sequential_g(form, ploughs)
#'  
#' out.boots <- boots_g(s1)
#' 
#' summary(out.boots)
#' 
#' @export
#'
boots_g <- function(seqg, boots = 1000) {
  stopifnot(boots %% 1 == 0 & boots > 0)
  
  B <- matrix(NA, 
              nrow = boots, 
              ncol = length(coef(seqg)),
              dimnames = list(NULL, names(coef(seqg))))
  nobs <- nrow(seqg$model)
  
  # init
  for (i in 1:boots) {
    if (i == 1) { 

    }
    
    # bootstrap model frame
    draw <- sample(1:nobs, replace = TRUE) 
    
    # run seqg
    seqg.draw <- sequential_g(formula = seqg$formula, 
                              data = seqg$model[draw, ])
     
    # store coef
    B[i, ] <- coef(seqg.draw)
  }
  
  class(B) <- "seqgboots"
  return(B)
}
