library(tidyverse)

AntColony <- function(f, city_coords, alpha = 1, beta = 4, rho = 0.5, m, n_iter) {
  
  n <- nrow(city_coords)
  
  if(n > m) {
    warning("m (the number of ants) are advised to be equal to, or bigger than n (the number of cities)")
  }
  
  
  ## Initialization
  # ant_route : a matrix whose row indicates ant in an iteration, and a column indicates a city
  # distance : a vector to store the distance an ant takes while walking all the cities.
  # adjacency : to be used for tau_increament.
  ant_route <- matrix(rep(0, n * m), nrow = m)
  distance <- rep(0, m)
  tau_increment <- adjacency <- weight <- distance <- matrix(rep(0, n * n), nrow = n)
  
  # distance : a matrix whose (i,j) th element indicates distance between i and j
  for(i in 1:n) {
    for(j in 1:n) {
      distance[i,j] <- sqrt((city_coords[i,1] - city_coords[j,1])^2 +
                           (city_coords[i,2] - city_coords[j,2])^2)
    }
  }
  
  # eta : a matrix whose (i,j) th element indicates approximity, 
  #       which is inversely proportaional to distance
  # Q : a constant to define eta; eta[i,j] = Q / distance.
  # We aim to set Q s.t. eta[i,j] ranges below 1.
  Q <- median(distance[row(distance) != col(distance)])
  eta <- Q / distance
  
  tau <- matrix(rep(1/n, n * n), nrow = n)
  
  
  ## iteration
  tau_list <- list()
  weight_list <- list()
  ant_route_opt_list <- list()
  distance_opt_list <- list()
  
  for(l in 1:n_iter) {
    # weight update
    for(i in 1:n) {
      for(j in 1:n) {
        weight[i,j] <- tau[i,j] ^ alpha * eta[i,j] ^ beta
      }
    }
    
    tour_length <- c()
    for(i in 1:m) {
      candidates <- 1:n
      for(j in 1:(n-1)) {
        prob <- c()
        if(j == 1) {
          ant_route[i,j] <- sample(1:n, 1)
        } else {
          prob_denom <- 0
          for(k in candidates) {
            prob_denom <- prob_denom + weight[ant_route[i,j-1],k]
          }
          for(k in seq_along(candidates)) {
            prob[k] <- weight[ant_route[i,j-1],candidates[k]] / prob_denom
          }
          ant_route[i,j] <- sample(candidates, 1, prob = prob)
        }
        candidates <- candidates %>% setdiff(ant_route[i,j])
      }
      ant_route[i,n] <- candidates
      
      tour_length[i] <- ant_route[i,] %>% as.vector() %>% f()
      
      # Calculating tau_increment for i-th iteration
      for(j in 1:n) {
        if(j == 1) {
          row <- ant_route[i,n] 
          col <- ant_route[i,j]
        } else {
          row <- ant_route[i,j-1]
          col <- ant_route[i,j]
        }
        adjacency[row,col] <- adjacency[row,col] + 1
      }
      tau_increment  <- tau_increment + adjacency / tour_length[i]
    }
    
    # tau update
    tau <- (1- rho) * tau + tau_increment 
    
    distance_opt <- min(tour_length)
    opt_ant_index <- which.min(tour_length)
    ant_route_opt <- ant_route[opt_ant_index,]
    
    # storing & before next iteration
    tau_list <- tau_list %>% append(list(tau))
    weight_list <- weight_list %>% append(list(weight))
    ant_route_opt_list <- ant_route_opt_list %>% append(list(ant_route_opt))
    distance_opt_list <- distance_opt_list %>% append(list(distance_opt))
    
    #clearing the adjacency and tau_increment matrices
    adjacency <- matrix(rep(0, n * n), nrow = n)
    tau_increment <- matrix(rep(0, n * n), nrow = n)
  }

  return(list(
    taus = tau_list, 
    weights = weight_list, 
    optimal_routes = ant_route_opt_list, 
    optimal_distances = distance_opt_list
  ))

}


### Example Data : 50 cities connection with Cartesian distance ----------
set.seed(123456)
cities <- data.frame(x = runif(50, max = 100), y = runif(50, max = 250))


### R function that return the total distance ----------
distance <- function(city_coords = cities, permutation) {
  dist <- 0
  n <- nrow(city_coords)
  city_coords <- city_coords[permutation, ]
  
  for(i in 1:(n-1)) {
    dist <- dist + sqrt((city_coords[i+1,1] - city_coords[i,1])^2 + 
                          (city_coords[i+1,2] - city_coords[i,2])^2)
  }
  dist <- dist + sqrt((city_coords[n,1] - city_coords[1,1])^2 + 
                        (city_coords[n,2] - city_coords[1,2])^2)
  
  return(dist)
}


### R function that return the total distance ----------
distance <- function(city_coords = cities, permutation) {
  dist <- 0
  n <- nrow(city_coords)
  city_coords <- city_coords[permutation, ]
  
  for(i in 1:(n-1)) {
    dist <- dist + sqrt((city_coords[i+1,1] - city_coords[i,1])^2 + 
                          (city_coords[i+1,2] - city_coords[i,2])^2)
  }
  dist <- dist + sqrt((city_coords[n,1] - city_coords[1,1])^2 + 
                        (city_coords[n,2] - city_coords[1,2])^2)
  
  return(dist)
}


### Rcpp function that return the total distance  ----------
library("Rcpp")

cppFunction('double distance_Rcpp(DataFrame city_coords, IntegerVector permutation) {
  double dist = 0;
  int nrow = city_coords.nrows();
  int ncol = 2;
  
  NumericVector v1 (nrow);
  NumericVector v2 (nrow);
  DataFrame permuted_coords = DataFrame::create(Named("v1") = v1, 
                                                Named("v2") = v2);
  
  for(int i = 0; i < ncol; ++i) {
    NumericVector col = city_coords[i];
    NumericVector new_col (nrow);
    
    for(int j = 0; j < nrow; ++j) {
      int index = permutation[j];
      new_col[j] = col[index-1];
    }
    permuted_coords[i] = new_col;
  }
  
  NumericVector permuted_x = permuted_coords[0];
  NumericVector permuted_y = permuted_coords[1];
    
  for(int i = 0; i < (nrow - 1); ++i) {
    dist = dist + sqrt(pow(permuted_x[i+1] - permuted_x[i], 2) + 
                       pow(permuted_y[i+1] - permuted_y[i], 2));
  }

  dist = dist + sqrt(pow(permuted_x[nrow-1] - permuted_x[0], 2) + 
                     pow(permuted_y[nrow-1] - permuted_y[0], 2));
  
  return(dist);
}')


### result  ----------
result <- AntColony(f = function(x) {distance_Rcpp(city_coords = cities, permutation = x)}, 
                    city_coords = cities, 
                    m = 1500, 
                    n_iter = 10)

### Plot : Optimal ----------

min_iter_index <-result$optimal_distances %>% which.min
opt_order <- result$optimal_routes %>% `[[`(min_iter_index)
route_optimal <- cities[opt_order, ]

route_optimal %>% ggplot(mapping = aes(x = x, y =y)) +
  geom_point() +
  geom_path()

### Animation : Route changes by Iteration----------
library(gganimate)

route_optimal <- tibble(
  x = numeric(),
  y = numeric(),
  iteration = integer()
)

for(i in 1:(result$optimal_routes %>% length())) {
  Run_order <- result$optimal_routes %>% `[[`(i)
  route_optimal_temp <- cities[Run_order, ]
  rounte_optimal_temp <- route_optimal_temp %>% bind_cols(iteration = i)
  route_optimal <- route_optimal %>% bind_rows(rounte_optimal_temp)
}

route_transition <- route_optimal %>% ggplot(mapping = aes(x = x, y = y)) +
  geom_point() +
  geom_path() +
  transition_states(iteration)

animate(route_transition)



### Session Info  ----------
library(sessioninfo)
session_info(include_base = TRUE)

# ─ Session info ────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31 ucrt)
# os       Windows 11 x64 (build 22621)
# system   x86_64, mingw32
# ui       RStudio
# language (EN)
# collate  Korean_Korea.utf8
# ctype    Korean_Korea.utf8
# tz       Asia/Seoul
# date     2024-04-05
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   NA

# ─ Packages ────────────────────────────────────────────────────────────────────────────────────
# ! package     * version  date (UTC) lib source
# base        * 4.3.2    2023-10-31 [?] local
# class         7.3-22   2023-05-03 [2] CRAN (R 4.3.2)
# classInt      0.4-10   2023-09-05 [1] CRAN (R 4.3.3)
# cli           3.6.2    2023-12-11 [1] CRAN (R 4.3.2)
# colorspace    2.1-0    2023-01-23 [1] CRAN (R 4.3.2)
# P compiler      4.3.2    2023-10-31 [2] local
# crayon        1.5.2    2022-09-29 [1] CRAN (R 4.3.2)
# P datasets    * 4.3.2    2023-10-31 [2] local
# DBI           1.2.2    2024-02-16 [1] CRAN (R 4.3.2)
# dplyr       * 1.1.4    2023-11-17 [1] CRAN (R 4.3.2)
# e1071         1.7-14   2023-12-06 [1] CRAN (R 4.3.3)
# fansi         1.0.6    2023-12-08 [1] CRAN (R 4.3.2)
# farver        2.1.1    2022-07-06 [1] CRAN (R 4.3.2)
# forcats     * 1.0.0    2023-01-29 [1] CRAN (R 4.3.2)
# generics      0.1.3    2022-07-05 [1] CRAN (R 4.3.2)
# gganimate   * 1.0.9    2024-02-27 [1] CRAN (R 4.3.3)
# ggplot2     * 3.5.0    2024-02-23 [1] CRAN (R 4.3.2)
# gifski        1.12.0-2 2023-08-12 [1] CRAN (R 4.3.3)
# glue          1.7.0    2024-01-09 [1] CRAN (R 4.3.2)
# P graphics    * 4.3.2    2023-10-31 [2] local
# P grDevices   * 4.3.2    2023-10-31 [2] local
# P grid          4.3.2    2023-10-31 [2] local
# gtable        0.3.4    2023-08-21 [1] CRAN (R 4.3.2)
# hms           1.1.3    2023-03-21 [1] CRAN (R 4.3.2)
# KernSmooth    2.23-22  2023-07-10 [2] CRAN (R 4.3.2)
# labeling      0.4.3    2023-08-29 [1] CRAN (R 4.3.1)
# lifecycle     1.0.4    2023-11-07 [1] CRAN (R 4.3.2)
# lpSolve       5.6.20   2023-12-10 [1] CRAN (R 4.3.2)
# lubridate   * 1.9.3    2023-09-27 [1] CRAN (R 4.3.2)
# magrittr      2.0.3    2022-03-30 [1] CRAN (R 4.3.2)
# P methods     * 4.3.2    2023-10-31 [2] local
# munsell       0.5.0    2018-06-12 [1] CRAN (R 4.3.2)
# pillar        1.9.0    2023-03-22 [1] CRAN (R 4.3.2)
# pkgconfig     2.0.3    2019-09-22 [1] CRAN (R 4.3.2)
# prettyunits   1.2.0    2023-09-24 [1] CRAN (R 4.3.2)
# progress      1.2.3    2023-12-06 [1] CRAN (R 4.3.2)
# proxy         0.4-27   2022-06-09 [1] CRAN (R 4.3.3)
# purrr       * 1.0.2    2023-08-10 [1] CRAN (R 4.3.2)
# R6            2.5.1    2021-08-19 [1] CRAN (R 4.3.2)
# Rcpp        * 1.0.12   2024-01-09 [1] CRAN (R 4.3.2)
# readr       * 2.1.5    2024-01-10 [1] CRAN (R 4.3.2)
# rlang         1.1.3    2024-01-10 [1] CRAN (R 4.3.2)
# rstudioapi    0.15.0   2023-07-07 [1] CRAN (R 4.3.2)
# scales        1.3.0    2023-11-28 [1] CRAN (R 4.3.3)
# sessioninfo * 1.2.2    2021-12-06 [1] CRAN (R 4.3.3)
# sf            1.0-16   2024-03-24 [1] CRAN (R 4.3.3)
# P stats       * 4.3.2    2023-10-31 [2] local
# stringi       1.8.3    2023-12-11 [1] CRAN (R 4.3.2)
# stringr     * 1.5.1    2023-11-14 [1] CRAN (R 4.3.2)
# tibble      * 3.2.1    2023-03-20 [1] CRAN (R 4.3.2)
# tidyr       * 1.3.1    2024-01-24 [1] CRAN (R 4.3.2)
# tidyselect    1.2.0    2022-10-10 [1] CRAN (R 4.3.2)
# tidyverse   * 2.0.0    2023-02-22 [1] CRAN (R 4.3.2)
# timechange    0.3.0    2024-01-18 [1] CRAN (R 4.3.2)
# P tools         4.3.2    2023-10-31 [2] local
# transformr    0.1.5    2024-02-26 [1] CRAN (R 4.3.3)
# tweenr        2.0.3    2024-02-26 [1] CRAN (R 4.3.3)
# tzdb          0.4.0    2023-05-12 [1] CRAN (R 4.3.2)
# units         0.8-5    2023-11-28 [1] CRAN (R 4.3.3)
# utf8          1.2.4    2023-10-22 [1] CRAN (R 4.3.2)
# P utils       * 4.3.2    2023-10-31 [2] local
# vctrs         0.6.5    2023-12-01 [1] CRAN (R 4.3.2)
# withr         3.0.0    2024-01-16 [1] CRAN (R 4.3.2)

# [1] C:/Users/GIHUN/AppData/Local/R/win-library/4.3
# [2] C:/Program Files/R/R-4.3.2/library

# P ── Loaded and on-disk path mismatch.

# ───────────────────────────────────────────────────────────────────────────────────────────────
# [1] C:/Users/GIHUN/AppData/Local/R/win-library/4.3
# [2] C:/Program Files/R/R-4.3.2/library

# P ── Loaded and on-disk path mismatch.
