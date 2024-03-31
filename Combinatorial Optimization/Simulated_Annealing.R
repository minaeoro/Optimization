library(tidyverse)
library(bench)
library(scales)

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
    
    # new solution
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
## Example

set.seed(123456)
cities <- data.frame(x = runif(50, max = 100), y = runif(50, max = 250))

# R function that return the total distance
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


# Rcpp function that return the total distance
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

# Benchmark : from package bench
mark_SA <- mark(R_function = SimulatedAnnealing(function(x) {distance(city_coords = cities, permutation = x)}, 50, max_iter = 10),
                Rccp_function = SimulatedAnnealing(function(x) {distance_Rcpp(city_coords = cities, permutation = x)}, 50, max_iter = 10),
                check = FALSE)

# Running Simulated Annealing : Run1 - by exponential; Run2 - by Constant Exponential 
Run1 <- SimulatedAnnealing(function(x) {distance_Rcpp(city_coords = cities, permutation = x)}, 50, max_iter = 200000, method = "exponential")
Run2 <- SimulatedAnnealing(function(x) {distance_Rcpp(city_coords = cities, permutation = x)}, 50, max_iter = 200000, method = "const_exponential")




# Plot : Total Distance at Each Iteration
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

    

library(gganimate)

# Plot : Optimal 
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

# Animation : Route changes by Iteration
route_transition <- route_optimal %>% ggplot(mapping = aes(x = x, y = y)) +
                    geom_point() +
                    geom_path() +
                    transition_states(iteration)

animate(route_transition)
