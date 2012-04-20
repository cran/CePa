# dif: vector, gene symbol
# bk: vector, gene symbol
# path: igraph object or edge list
# mapping: map gene symbol to node id, matrix or data.frame
#          node.id  gene.id
# cen: centrality method
# cen.name: centrality name
cepa = function(dif, bk, pathway, mapping = cbind(bk, bk), cen = "equal.weight",
                cen.name = if(is.function(cen)) deparse(substitute(cen)) else if(mode(cen) == "name") deparse(cen) else cen,
                iter = 1000) {

    if(length(dif) > length(bk)) {
        stop("Length of differential genes should not be larger than the length of background genes.\n")
    }
    if(sum(dif %in% bk) != length(dif)) {
        stop("Differential genes must be all in background list.\n")
    }
    if(is.matrix(pathway) || is.data.frame(pathway)) {
        if(length(dim(pathway)) != 2 || dim(pathway)[2] != 2) {
            stop("if pathway is a matrix or data frame, it should be 2 dimension and the number of columns is 2.\n")
        }
        pathway = generate.pathway(pathway)
    }
    if(sum(dif %in% mapping[, 2]) == 0) {
        stop("Cannot map gene to node. Check the mapping argument.\n")
    }
    if(iter < 100) {
        stop("Iterations should not be smaller than 100.\n")
    }
    
    if(length(cen) > 1) {
        stop("Length of cen must be equal to 1.\n") 
    }
    
    weight = centrality(pathway, cen)
    
    if(any(weight < 0)) {
        stop("Weight should not be negative.")
    }

    add = 0
    # if there are none-zero weight values
    if(any(weight > 0)) {
        add = min(weight[weight > 0])/100
    }
    weight = weight + add
    
    # nodes in the pathway
    node = pathway.nodes(pathway)
    
    # only the mapping in the pathway
    mapping = mapping[mapping[, 1] %in% node, ]
    
    # get node names formatted with genes
    node.name = node
    node.diff.gene = character(length(node))
    member = character(0)
    for(i in 1:length(node)) {
        # genes that exsit in the node
        l = mapping[, 1] == node[i]
        
        # if find nodes with genes mapped
        if(sum(l)) {
            member = sort(unique(mapping[l, 2]))
            # mark the diff genes
            node.diff.gene[i] = paste(member[member %in% dif], collapse = ", ")
            member[member %in% dif] = paste("[", member[member %in% dif], "]", sep="")
            node.name[i] = paste(member, collapse = "\n")
        }
    }
    
    # map dif genes to node id
    dif.node = unique(mapping[mapping[, 2] %in% dif, 1])
    
    is.dif.node = as.numeric(node %in% dif.node)
    s = sum(is.dif.node * weight)
    if(sum(is.dif.node > 0) == 0)
        ds = c(0, 0, 0, 0)
    else
        ds = quantile(weight[is.dif.node > 0], c(1, 0.75, 0.5, 0))
    names(ds) = c("max", "q75", "median", "min")
    
    # sampling
    p.dif = length(dif) / length(bk)
    s.random = numeric(iter)
    ds.random = matrix(0, iter, 4)   # descriptive statistic of the node
    
    # genes in the pathway
    gene = unique(mapping[mapping[, 1] %in% node, 2])
    
    colnames(ds.random) = c("max", "q75", "median", "min")
    for(i in 1:iter) {
        # first select from pathway genes
        dif.random = gene[as.logical(rbinom(length(gene), 1, p.dif))]
        
        # then map to node id
        dif.node.random = unique(mapping[mapping[, 2] %in% dif.random, 1])
        # find which node is differentially affected
        is.dif.node.random = as.numeric(node %in% dif.node.random)
        # calculate the score
        s.random[i] = sum(is.dif.node.random * weight)
        if(sum(is.dif.node.random > 0) == 0) {
            ds.random[i, ] = c(0, 0, 0, 0)
        }
        else {
            ds.random[i, ] = quantile(weight[is.dif.node.random > 0], c(1, 0.75, 0.5, 0))
        }
    }
    
    p.value = (sum(s.random >= s) + 1) / (iter + 1)
    
    # calculate fisher's exact test
    # 2x2 contingency table for fisher's exact test
    dif.gene = intersect(dif, gene)
    m = matrix(0, 3, 3)
    m[1, 1] = length(dif.gene)
    m[3, 1] = length(gene)
    m[2, 1] = m[3, 1] - m[1, 1]
    m[1, 3] = length(dif)
    m[1, 2] = m[1, 3] - m[1, 1]
    m[3, 3] = length(bk)
    m[2, 3] = m[3, 3] - m[1, 3]
    m[3, 2] = m[3, 3] - m[3, 1]
    m[2, 2] = m[3, 2] - m[1, 2]
    # one-side test
    p.ora = fisher.test(m[1:2, 1:2], alternative="greater")$p.value
    
    n.dif.node = length(dif.node)
    n.dif.gene = length(dif.gene)
    n.node = length(node)
    n.gene = length(gene)
    
    count = c(n.dif.node, n.node, n.dif.gene, n.gene)
    names(count) = c("n.dif.node", "n.node", "n.dif.gene", "n.gene")
    
    res = list("score" = s,
               "ds" = ds,
               "ds.simulation" = ds.random,
               "p.value" = p.value,
               "p.ora" = p.ora,
               "simulation" = s.random,
               "centrality" = cen.name,
               "weight" = weight,
               "is.dif.node" = as.logical(is.dif.node),
               "node.name" = node.name,
               "node.diff.gene" = node.diff.gene,
               "count" = count,
               "pathway" = pathway)
    
    class(res) = "cepa"
    
    return(invisible(res))

}



