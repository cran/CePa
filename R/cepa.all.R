
# calculate all pathway's significance with all centrality provided
# dif: diff genes   gene_id
# bk: background genes     gene_id
# pathList: list, should have names
# interactionList: three column matrix
#                  interaction_id   input_id    output_id
# mapping: two column matrix
#          node_id    gene_id
# centrality list: if use functions, use quote or substitute
# iter: times of simulations
cepa.all = function(dif, bk, pathList, interactionList, mapping = cbind(bk, bk),
           cen = c("equal.weight", "in.degree", "out.degree", "betweenness", "in.reach", "out.reach"),
           cen.name = sapply(cen, function(x){if(mode(x) == "name") deparse(x) else x}), iter = 1000,
           min.node = 5, max.node = 500) {
    
    if(!is.list(pathList)) {
        stop("pathList should be a list.\n")
    }
    if(is.null(names(pathList))) {
        stop("pathList should have names.\n")
    }
    if(length(dim(interactionList)) != 2) {
        stop("interactionList should be two dimension matrix.\n")
    }
    if(dim(interactionList)[2] != 3) {
        stop("interactinList should contain 2 columns.\n")
    }
    if(length(cen) < 1) {
        stop("cen argument must be specified.\n")
    }
    for(ce in cen) {
        if(is.function(ce)) {
            stop("Functions cannot be used directly, use quote or substitute.\n")
        }
    }
    
    l = sapply(pathList, function(x) {
                   it = interactionList[interactionList[, 1] %in% x, 2:3]
                   node = unique(c(it[, 1], it[, 2]))
				   l.node = length(node)
				   gene = unique(mapping[mapping[,1] %in% node, 2])
				   l.gene = length(gene)
				   return(c(l.node, l.gene))
               })
    pathList = pathList[l[1, ] >= min.node & l[1, ] <= max.node
	                    | l[2, ] >= min.node & l[2, ] <= max.node]
    
    pathway.name = names(pathList)    
    pathway.result = list()
    length(pathway.result) = length(cen)
    pathway.result = lapply(pathway.result, function(x) {
                              y = list()
                              length(y) = length(pathway.name)
                              names(y) = pathway.name
                              return(y)
                            })
    names(pathway.result) = cen.name
    
    for(i in 1:length(pathList)) {
        
        cat("  ", i, "/", length(pathList), ", ", pathway.name[i], "...\n", sep="")
        
        # interaction id in pathway
        path = pathList[[i]]
        inter = interactionList[interactionList[, 1] %in% path, 2:3]
        # generate graph from edge list
        pathway = generate.pathway(as.matrix(inter))
        
        j = 0
        for(ce in cen) {
            j = j + 1
            pathway.result[[j]][[i]] = cepa(dif, bk, pathway, mapping, ce, iter=iter)
        }
    }

    res = list("pathway.name" = pathway.name,
               "pathway.result" = pathway.result)
    
    class(res) = "cepa.all"
    return(res)
}
