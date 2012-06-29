# p values of all pathways
p.table = function(x) {

    if(class(x) != "cepa.all") {
        stop("x should be cepa.all object.\n")
    }
    
    n.pathway = length(x)
    p.value = matrix(0, nrow=length(x), ncol= length(x[[1]]))
    for(i in 1:length(x)) {
        p.value[i, ] = sapply(x[[i]], function(x) x$p.value)
    }

    rownames(p.value) = names(x)
    colnames(p.value) = names(x[[1]])
    
    return(p.value)
}
