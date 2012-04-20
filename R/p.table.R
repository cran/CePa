
# p values of all pathways
p.table = function(x) {

    if(class(x) != "cepa.all") {
        stop("x should be cepa.all object.\n")
    }
    
    n.centrality = length(x$pathway.result)
    p.value = matrix(0, nrow=length(x$pathway.name), ncol= n.centrality)
    for(i in 1:n.centrality) {
        p.value[, i] = sapply(x$pathway.result[[i]], function(x) x$p.value)
    }

    rownames(p.value) = x$pathway.name
    colnames(p.value) = names(x$pathway.result)
    
    return(p.value)
}
