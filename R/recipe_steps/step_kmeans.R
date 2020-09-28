

step_kmeans <- function(
  recipe, ..., 
  role = NA, 
  trained = FALSE, 
  centers = 5,
  objects = NULL,
  options = list(),
  skip = FALSE,
  id = rand_id("kmeans")
) {

  add_step(
    recipe, 
    step_kmeans_new(
      terms = terms <- ellipse_check(...), 
      trained = trained,
      role = role, 
      centers = centers,
      objects = objects,
      options = options,
      skip = skip,
      id = id
    )
  )
}


step_kmeans_new <- 
  function(terms, role, trained, centers, objects, options, skip, id) {
    step(
      subclass = "kmeans", 
      terms = terms,
      role = role,
      trained = trained,
      centers = centers,
      objects = objects,
      options = options,
      skip = skip,
      id = id
    )
  }


prep.step_kmeans <- function(x, training, info = NULL, ...){
  col_names <- terms_select(terms = x$terms, info = info) 
 
  objects <- purrr::map(training[, col_names], kmeans, centers = x$centers)
  
  step_kmeans_new(
    terms = x$terms, 
    trained = TRUE,
    role = x$role, 
    centers = x$centers,
    objects = objects,
    options = x$options,
    skip = x$skip,
    id = x$id
  )
}


assign_cluster <- function(x, centers) {
  # compute squared euclidean distance from each sample to each cluster center
  tmp <- sapply(seq_len(length(x)),
                function(i) apply(centers, 1,
                                  function(v) sum((x[i]-v)^2)))
  factor(max.col(-t(tmp)), levels = 1:nrow(centers))  # find index of min distance
}


bake.step_kmeans <- function(object, new_data, ...) {
  require(tibble)
  
  for (i in names(object$objects)) {
    centers <- object$objects[[i]]$centers
    new_data[, i] <- assign_cluster(getElement(new_data, i), 
                                    centers[order(centers[, 1]), , drop = FALSE])
  }
    
  as_tibble(new_data)
}


print.step_kmeans <- function (x, width = max(20, options()$width - 35), ...) 
{
  cat("K-means discretization with", x$centers, "centers")
  if (x$trained) {
    cat(" [trained]\n")
  }
  else {
    cat("\n")
  }
  invisible(x)
}




