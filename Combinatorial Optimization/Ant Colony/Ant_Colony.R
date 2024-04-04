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