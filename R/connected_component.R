is_connected <- function(fw){
  
  DSF = function(m, i, visited){
    if (visited[i]) {
      return(1)
    }else{
      visited[i] = TRUE
      for (n in which(m[i,] >0)){
        DSF(m, n, visited)
      }
      }
  }
  
  # make the network undirected to simply use DSF algorithm
  fw.s <- fw
  fw.s[lower.tri(fw.s)] = t(fw.s)[lower.tri(fw.s)]
  # create the boolean vector of visited nodes
  visited = rep(FALSE, nrow(fw.s))
  # call deep search first alg. on the first node (should it be randomly chosen?)
  DSF(fw.s, 1, visited)
  
  # check if all nodes where reach by the DSF, and return:
  return(all(visited))
}