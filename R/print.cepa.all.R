# summary for top.all object
print.cepa.all = function(x, ...) {
    cat("\n")
    cat("number of pathways:", length(x$pathway.name), "\n")
    cat("\n")
    
    p.value = p.table(x)
    centrality = colnames(p.value)

    cat("Significant pathways (p.value <= 0.01):\n")
    sig = matrix(0, nrow=length(centrality), 1)
    rownames(sig) = centrality
    colnames(sig) = "Number"
    for(i in 1:length(centrality)) {
        sig[i, 1] = sum(p.value[ ,i] <= 0.01)
    }
    print(sig)
    
    cat("\n")
}
