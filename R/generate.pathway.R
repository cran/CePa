# construct the igraph object pathway
generate.pathway = function(el) {
    
    if(dim(el)[2] != 2) {
        stop("Second dimension of edgelist should be 2.\n")
    }
    
    el = as.matrix(el)
    
    g = graph.edgelist(el, directed = TRUE)
    
    return(g)
}
