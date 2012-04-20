print.cepa = function(x, ...) {
    cat("\n")
    cat("  diff nodes in pathway:", x$count["n.dif.node"], "\n")
    cat("  all nodes in pathway:", x$count["n.node"], "\n")
    cat("  diff genes in pathway:", x$count["n.dif.gene"], "\n")
    cat("  all genes in pathway:", x$count["n.gene"], "\n")
    cat("\n\n")
    cat("  weight:", x$centrality, "\n")
    cat("  p-value:", sprintf("%.1e", x$p.value), "\n")
    cat("\n")
}