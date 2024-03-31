library(tidyverse)

### Simulated Annealing Function ----------
# x : a vector, which is a permutation of 1:n
# f : function that takes x and return.
# Objective : finding the minimum
SimulatedAnnealing <- function(f, n, max_iter, method = "exponential") {
  x <- 1:n  #initial solution 
  y <- f(x) 
  x_best <- x; y_best <- y
  iter = c() # n-th iteration
  x_optim = c()
  y_optim = c() 
  y_current = c()
  
  temp <- 1000 # Initial Temperature
  rate <- exp(log(0.1/temp)/max_iter) # Exponential scheme : chosen s.t. after iterations of max_iter, the temperature becomes 0.1
  
  for(k in 1:max_iter) {
    
    # new solution : two points are selected, after which all the points in-between (along with the two points) are reversed.
    pivots <- sample(1:n, 2, replace = FALSE)
    pivot_smaller <- min(pivots)
    pivot_bigger <- max(pivots)
    
    x_rev_part <- rev(x[pivot_smaller:pivot_bigger])
    
    if(pivot_smaller > 1 & pivot_bigger < n) {
      x_non_rev_front_part <- x[1:(pivot_smaller-1)]
      x_non_rev_back_part <- x[(pivot_bigger+1):n]
    } else if(pivot_smaller == 1 & pivot_bigger < n) {
      x_non_rev_front_part <- c()
      x_non_rev_back_part <- x[(pivot_bigger+1):n]
    } else if(pivot_smaller > 1 & pivot_bigger == n) {
      x_non_rev_front_part <- x[1:(pivot_smaller-1)]
      x_non_rev_back_part <- c()
    } else {
      x_non_rev_front_part <- c()
      x_non_rev_back_part <- c()
    }
    
    x_new <- c(x_non_rev_front_part, x_rev_part, x_non_rev_back_part)
    
    # new distance
    y_new <- f(x_new)
    
    y_new_diff <- y_new - y
    
    if(y_new_diff < 0 | runif(1) < exp(-y_new_diff/temp)) {
      x <- x_new
      y <- y_new 
    }
    
    if(y_new < y_best) {
      x_best <- x_new
      y_best <- y_new
    }
    
    iter = c(iter, k)
    x_optim = c(x_optim, paste(x_best, collapse = " "))
    y_optim = c(y_optim, y_best)
    y_current = c(y_current, y_new)
    
    # Temperature Scheme : Exponential or Constant Exponential ?
    if(method == "exponential") {
      temp <- rate * temp
    } else if(method == "const_exponential") {
      if(k < max_iter / 5) {
        next
      } else {
        temp <- rate * temp
      }
    }
    else {
      stop("Wrong method name : it should be either 'exponential' or 'const_exponential.'")
    }
  }
  
  return(list(x_opt = x_best,
              y_opt = y_best,
              records = tibble(iter = iter,
                               x_optim = x_optim,
                               y_optim = y_optim,
                               y_current = y_current)))
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

### Benchmark : from package bench ----------
library(bench)
mark_SA <- mark(R_function = SimulatedAnnealing(function(x) {distance(city_coords = cities, permutation = x)}, 50, max_iter = 10),
                Rccp_function = SimulatedAnnealing(function(x) {distance_Rcpp(city_coords = cities, permutation = x)}, 50, max_iter = 10),
                check = FALSE)

mark_SA %>% plot(type = "violin")

# Running Simulated Annealing : Run1 - by exponential; Run2 - by Constant Exponential 
Run1 <- SimulatedAnnealing(function(x) {distance_Rcpp(city_coords = cities, permutation = x)}, 50, max_iter = 200000, method = "exponential")
Run2 <- SimulatedAnnealing(function(x) {distance_Rcpp(city_coords = cities, permutation = x)}, 50, max_iter = 200000, method = "const_exponential")




### Plot : Total Distance at Each Iteration ----------
library(scales)
y_records <- Run2$records %>% select(iter, y_optim, y_current)

y_records %>% 
  ggplot(mapping = aes(x = iter)) +
    geom_line(mapping = aes(y = y_optim), color = "tomato", linewidth = 1.2) +
    geom_line(mapping = aes(y = y_current), color ="black") +
    theme_minimal() +
    labs(x = "Iterations",
         y = "Total Distance",
         color = "Value",
         title = "Total Distance at Each Iteration",
         caption = "Red Line : Minimum Distance at the iteration.\nBlack Line : Current Distance found at the iteration.") +
    scale_x_continuous(labels = scales::comma) +  
    scale_y_continuous(labels = scales::comma)


### Plot : Optimal ----------
Run_order <- Run2$records$x_optim %>% tail(1) %>% strsplit(" ") %>% unlist() %>% as.numeric()
route_optimal <- cities[Run_order, ]

route_optimal %>% ggplot(mapping = aes(x = x, y =y)) +
  geom_point() +
    geom_path()


route_optimal <- tibble(
  x = numeric(),
  y = numeric(),
  iteration = integer()
)

for(i in c(1, 500 * (1:400))) {
  Run_order <- Run2$records$x_optim[i] %>% strsplit(" ") %>% unlist() %>% as.numeric()
  route_optimal_temp <- cities[Run_order, ]
  rounte_optimal_temp <- route_optimal_temp %>% bind_cols(iteration = i)
  route_optimal <- route_optimal %>% bind_rows(rounte_optimal_temp)
}

### Animation : Route changes by Iteration----------
library(gganimate)
route_transition <- route_optimal %>% ggplot(mapping = aes(x = x, y = y)) +
                    geom_point() +
                    geom_path() +
                    transition_states(iteration)

animate(route_transition)


### Session Info  ----------
library(sessioninfo)
session_info(include_base = TRUE)


# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31 ucrt)
# os       Windows 11 x64 (build 22621)
# system   x86_64, mingw32
# ui       RStudio
# language (EN)
# collate  Korean_Korea.utf8
# ctype    Korean_Korea.utf8
# tz       Asia/Seoul
# date     2024-04-01
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   NA

# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ! package     * version  date (UTC) lib source
# base        * 4.3.2    2023-10-31 [?] local
# bench       * 1.1.3    2023-05-04 [1] CRAN (R 4.3.3)
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
# profmem       0.6.0    2020-12-13 [1] CRAN (R 4.3.3)
# progress      1.2.3    2023-12-06 [1] CRAN (R 4.3.2)
# proxy         0.4-27   2022-06-09 [1] CRAN (R 4.3.3)
# purrr       * 1.0.2    2023-08-10 [1] CRAN (R 4.3.2)
# R6            2.5.1    2021-08-19 [1] CRAN (R 4.3.2)
# Rcpp        * 1.0.12   2024-01-09 [1] CRAN (R 4.3.2)
# readr       * 2.1.5    2024-01-10 [1] CRAN (R 4.3.2)
# rlang         1.1.3    2024-01-10 [1] CRAN (R 4.3.2)
# rstudioapi    0.15.0   2023-07-07 [1] CRAN (R 4.3.2)
# scales      * 1.3.0    2023-11-28 [1] CRAN (R 4.3.3)
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
