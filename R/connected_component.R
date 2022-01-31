#' @title Detect whether a food web is composed of several disconected sub-networks
#'
#' @description Run a deep search first algorithm (DFS)
#'
#' @param fw binary adjacency matrix of the food web.
#'
#' @return Boolean: TRUE if the food web is connected, FALSE if several disconnected sub-networks are detected.
#'
#' @examples
#' 
#' library(ATNr)
#' set.seed(123)
#' # number of species, nutrients, and body masses
#' n_species <- 20
#' n_basal <- 5
#' n_nutrients <- 3
#' masses <- sort(10^runif(n_species, 2, 6)) #body mass of species
#' # create food web matrix
#' L <- create_Lmatrix(masses, n_basal)
#' L[, 1:n_basal] <- 0
#' fw <- L
#' fw[fw > 0] <- 1
#' connected = is_connected(fw)

is_connected <- function(fw){
  
  DFS = function(m, i, visited, env = envf){
    # loop over successors
    for (n in which(m[i,] > 0)){
      if (!env$visited[n]){
        # if n not visited, apply recursivity
        env$visited[n] = TRUE
        DFS(m, n, env$visited)
      } 
    }
    return()
  }
  
  # make the network undirected to simply use DSF algorithm
  fw.s <- fw
  fw.s[lower.tri(fw.s)] = t(fw.s)[lower.tri(fw.s)]
  # create the Boolean vector of visited nodes
  visited = rep(FALSE, nrow(fw.s))
  visited[1] = TRUE
  # call deep first search alg. on the first node (should it be randomly chosen?)
  envf = environment()
  visited <- DFS(fw.s, 1, visited, envf)
  
  # check if all nodes where reach by the DSF, and return:
  return(all(visited))
}
