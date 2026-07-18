### LIBRARIES ###
library(readxl)
library(tidyverse)

### FUNCTION DEV ###
#' zprime: calculates the Z-like directional statistic Z-prime using
#' Gamma-Poisson Conjugacy
#'
#' Follows derivation seen in zprime_derivation.pdf
#'
#' @param tada_data gene-level data including:
#' counts, mutation rates, and hyperparameters 
#' @param sheet_num excel sheet number
#' 
#' @return data-frame of Z-prime values for all genes 
#' 
#' @examples zprime("satterstrom_counts.xlsx", 2)
zprime <- function(tada_data, sheet_num){
  
  # read in data and extract relevant columns
  data <- read_excel(paste0("../data/raw/", tada_data), sheet = sheet_num, col_names = TRUE)
  
  input <- data %>%
    dplyr::select(gene, dn.ptv, mut.ptv, gamma_dn.ptv) %>%
    dplyr::rename("x" = "dn.ptv",
                  "mu" = "mut.ptv",
                  "gamma_bar" = "gamma_dn.ptv")
  
  # set constants
  beta <- 0.2 # dispersion/rate hyperparameter
  Pi <- 0.05 # ASD risk proportion
  n <- 6430 # number of trios
  gamma_bar <- input$gamma_bar # smoothed relative risk hyperparameter

  
  # relative risk coefficient
  c_i <- 2*n*(input$mu)
  
  # counts 
  x_i <- input$x
  
  # null marginal likelihood
  null_marg <- dpois(input$x, c_i)
  
  # alternative marginal likelihood
  alt_marg <- dnbinom(x_i, gamma_bar*beta, beta/(beta + c_i))
  
  # weights
  w_i <- Pi*alt_marg/(Pi*alt_marg + (1-Pi)*null_marg)
  
  # weighted posterior mean
  weighted_mean <- (1 - w_i) + w_i*((gamma_bar*beta + x_i)/(beta + c_i))
  
  # weighted posterior variance
  weighted_variance <- w_i*((gamma_bar*beta + x_i)/(beta + c_i)^2) + w_i*(1-w_i)*(((gamma_bar*beta + x_i)/(beta + c_i)) - 1)^2
  
  # final z-prime
  z_prime <- (weighted_mean - 1)/sqrt(weighted_variance)
  
  # get gene names ane z-prime
  result <- input %>% 
    cbind(z_prime) %>%
    dplyr::select(gene, z_prime)
  

  return(result)
}
