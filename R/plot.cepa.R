# node.attributes
plot.cepa = function(x, node.name = NULL, node.type = NULL,
                     graph.node.max.size = 20, graph.node.min.size = 3, graph.layout.method = NULL, ...) {
    
    s = x$score
    s.random = x$simulation
    ds = x$ds
    ds.random = x$ds.simulation
    weight = x$weight
    centrality = x$centrality
    p.value = x$p.value
    is.dif.node = x$is.dif.node
    pathway = x$pathway              # igraph object
    
    n.node = x$count["n.node"]
    
    node.id = pathway.nodes(pathway)
    # check node.name and node.type
    if(! is.null(node.name)) {
        node.name = node.name[node.id]
        if(sum(is.na(node.name))) {
            stop("cannot match all nodes with node name.\n")
        }
    }
    if(! is.null(node.type)) {
        node.type = node.type[node.id]
        if(sum(is.na(node.type))) {
            stop("cannot match all nodes with node type.\n")
        }
        if(length(setdiff(node.type, c("protein", "complex", "family", "compound", "rna")))) {
            stop("node.type only permits protein, complex, family, subunit, compound, rna.\n")
        }
    }
    
    if(length(weight) == 0) {
        plot(1,1, axes=FALSE, ann=FALSE, type="n")
        text(1,1,"empty pathway!", cex=3)
        return(NULL)
    }

    # color settigs
    color = character(4)
    names(color) = colnames(ds.random)
    color["max"] = rgb(227, 26, 28, maxColorValue = 255)
    color["q75"] = rgb(255, 191, 111, maxColorValue = 255)
    color["median"] = rgb(127, 201, 127, maxColorValue = 255)
    color["min"] = rgb(31, 120, 180, maxColorValue = 255)
    
    layout(rbind(c(1, 2, 4),
                 c(3, 0, 4)),
            widths=c(4, lcm(3), 8))
    
    # figure A
    mar = c(5.1, 4.1, 4.1, 2.1)
    par(mar = mar+c(0, 0, 0, -2))
    xrange = range(c(s.random, s))
    yrange = range(weight)
    if(yrange[1] == yrange[[2]]) {
        yrange = sort(c(0, 2*yrange[1]))
    }
    matplot(jitt(s.random, xrange), jitt(ds.random, yrange), xlim = xrange, ylim = yrange, col=color, pch=20, cex=0.2, xlab = "Simulated score", ylab="Centrality", main = paste("(A) Distribution of", centrality, "centrality\nin the pathway under simulation"))
    points(rep(s, 4), jitt(ds, yrange), col=color, pch=20, cex=5)
    legend(min(xrange), max(yrange), colnames(ds.random), col=color, pch=20, pt.cex=2)
    abline(v = s, col="red", lwd = 2)
    par(mar=mar)
    
    # figure B
    mar = c(5.1, 4.1, 4.1, 2.1)
    par(mar = mar+c(0, -2, 0, 2))
    plot(c(0.3, 1.7), range(weight), axes=FALSE, ann=FALSE, type="n", xaxs="i")
    axis(4)
    points(runif(length(weight), min=0.5, max=1.5), jitt(weight), pch=20, cex=0.4)
    mtext("Centrality of all nodes", 4, line=2, cex=0.8)
    title("(B)")
    par(mar=mar)
    
    # figure C
    mar = c(5.1, 4.1, 4.1, 2.1)
    par(mar = mar+c(0, 0, 0, -2))
    hist(s.random, freq = FALSE, xlim = range(c(s.random, s)), breaks = 50, xlab = "Simulated score", main = paste("(C) Histogram of simulated scores in the pathway\nusing", centrality, "centrality as weight"))
    box()
    abline(v = s, col = "red", lwd = 2)
    par(mar=mar)
    
    # figure D, network grpah
    # node should only be in protein, complex, family, compound, rna, subunit
    if(is.null(node.type)) {
    
        v.color = rep(rgb(145, 207, 96, 128, maxColorValue = 255), n.node)
        v.color[is.dif.node] = rgb(215, 48, 39, 128, maxColorValue = 255)
        complex.id = grepl("^\\d+$", x$node.name)
        v.color[complex.id] = rgb(50, 136, 189, 128, maxColorValue = 255)
        
        # shape
        v.shape = rep("circle", length(weight))
        v.shape[complex.id] = "square"
        
        v.frame.color = "white"
    }
    else {    

        v.color = rep(rgb(145, 207, 96, 128, maxColorValue = 255), n.node)
        v.color[node.type == "protein"] = rgb(145, 207, 96, 128, maxColorValue = 255)
        v.color[node.type == "complex"] = rgb(145, 207, 96, 128, maxColorValue = 255)
        v.color[node.type == "family"] = rgb(145, 207, 96, 128, maxColorValue = 255)
        v.color[node.type == "subunit"] = rgb(145, 207, 96, 128, maxColorValue = 255)
        v.color[node.type == "compound"] = rgb(50, 136, 189, 128, maxColorValue = 255)
        v.color[node.type == "rna"] = rgb(84, 39, 136, 128, maxColorValue = 255)
        v.color[is.dif.node] = rgb(215, 48, 39, 128, maxColorValue = 255)
        
        v.shape = rep("circle", length(node.id))
        v.shape[node.type == "protein"] = "circle"
        v.shape[node.type == "complex"] = "circle"
        v.shape[node.type == "family"] = "circle"
        v.shape[node.type == "subunit"] = "circle"
        v.shape[node.type == "compound"] = "square"
        v.shape[node.type == "rna"] = "square"
        
        v.frame.color = rep("white", length(node.id))
        v.frame.color[node.type == "protein"] = "white"
        v.frame.color[node.type == "complex"] = rgb(228, 26, 28, maxColorValue = 255)
        v.frame.color[node.type == "family"] = rgb(55, 126, 184, maxColorValue = 255)
        v.frame.color[node.type == "subunit"] = rgb(152, 78, 163, maxColorValue = 255)
    }
    
    # size
    max.size = graph.node.max.size
    min.size = graph.node.min.size
    v.size = rep(min.size, n.node)
    if(max(weight) != min(weight)) {
        v.size = min.size + (max.size - min.size)/(max(weight) - min(weight))*(weight - min(weight))
    }
    
    # label
    if(is.null(node.name)) {
        v.label = x$node.name
    }
    else {
        v.label = node.name
        v.label = gsub("/", "\n", v.label)
    }
    
    V(pathway)$color = v.color
    V(pathway)$shape = v.shape
    V(pathway)$size = v.size
    V(pathway)$label = v.label
    V(pathway)$frame.color = v.frame.color
    V(pathway)$label.cex = 1
    V(pathway)$label.font = 2
    V(pathway)$label.color = "black"
    E(pathway)$arrow.size = 0.5
    E(pathway)$color= "black"
    
    V(pathway)$node.diff.gene = x$node.diff.gene
    
    mar = c(5.1, 4.1, 4.1, 2.1)
    par(mar = mar+c(2, 0, 0, 0))
    if(is.null(graph.layout.method)) {
        layout.method = layout.reingold.tilford
    }
    else {
        layout.method = graph.layout.method
    }
    
    plot.igraph2(pathway, layout.method=layout.method)
    
    if(is.null(node.type)) {
        legend(-0.08, 0, c("Diff gene/complex/family", "gene/complex/family", "other type of nodes"),
               col=c(rgb(215, 48, 39, 128, maxColorValue = 255),
                     rgb(145, 207, 96, 128, maxColorValue = 255),
                     rgb(50, 136, 189, 128, maxColorValue = 255)), pch=c(16,16,15), pt.cex=2, yjust=0)
    }
    else {
        legend(-0.08, 0, c("Differential nodes", "Non-diff nodes", "complex", "family", "subunit", "compound", "rna"),
               col=c(rgb(215, 48, 39, 128, maxColorValue = 255),
                     rgb(145, 207, 96, 128, maxColorValue = 255),
                     rgb(228, 26, 28, maxColorValue = 255),
                     rgb(55, 126, 184, maxColorValue = 255),
                     rgb(152, 78, 163, maxColorValue = 255),
                     rgb(50, 136, 189, 128, maxColorValue = 255),
                     rgb(84, 39, 136, 128, maxColorValue = 255)),
                pch=c(16,16,21,21,21,15,15), pt.cex=c(2.5,2.5,2,2,2,2.5,2.5), cex=1.2, yjust=0)
    }
    title("(D) Graph view of the pathway")
    par(mar=mar, xpd=FALSE)
    layout(matrix(1, 1, 1))
    
    return(invisible(pathway))
}


# add random noise to a point in a xy coordinate system
# similar to function jitter
jitt = function(x, r = range(x[x != Inf], na.rm=TRUE)) {
    return(x + runif(length(x), -0.5, 0.5)*(r[2] - r[1])/50)
}

plot.igraph2 = function(g, layout.method = layout.random, ...) {
	
	v.color = V(g)$color
	v.shape = V(g)$shape
	v.shape2 = rep(16, length(v.shape))
	v.shape2[v.shape == "square"] = 15
	v.shape = v.shape2
	v.size = V(g)$size
	names(v.size) = pathway.nodes(g)
	v.label = V(g)$label
	v.frame.color = V(g)$frame.color
	v.label.cex = V(g)$label.cex
	v.label.font = V(g)$label.font
	v.label.color = V(g)$label.color
	e.arrow.size = E(g)$arrow.size
	e.color = E(g)$color
	
	ly = layout.method(g)
	r1 = max(ly[, 1]) - min(ly[, 1])
	ly[, 1] = (ly[, 1] - min(ly[, 1]))/r1
	r2 = max(ly[, 2]) - min(ly[, 2])
	ly[, 2] = (ly[, 2] - min(ly[, 2]))/r2
	rownames(ly) = pathway.nodes(g)
	# points
	plot(ly, cex = v.size,
	         col = v.color,
			 pch = v.shape,
			 ann = FALSE,
			 axes = FALSE,
			 asp = 1,
			 xlim = c(-0.1, 1.1),
			 ylim = c(-0.1, 1.1),
			 xaxs = "i",
			 yaxs = "i"
			 )
	box()
	frame.pch = v.shape
	frame.pch[frame.pch == 16] = 21
	frame.pch[frame.pch == 15] = 0
	points(ly, pch = frame.pch, cex = v.size,
	                     col = v.frame.color)
	text(ly, v.label, cex = v.label.cex,
					   col = v.label.color)
                       					   
	# edges
	el = get.edgelist(g)
	from = ly[el[, 1], ]
	from.node = rownames(from)
	to = ly[el[, 2], ]
	to.node = rownames(to)
	new.from = matrix(0, nrow=dim(from)[1], ncol=2)
	new.to = matrix(0, nrow=dim(to)[1], ncol=2)
	for(i in 1:length(from.node)) {
		foo = reedge(from[i, 1], from[i, 2], to[i, 1], to[i, 2], v.size[from.node[i]]/36/par("pin")[1], v.size[to.node[i]]/36/par("pin")[1])
		new.from[i, ] = foo[1:2]
		new.to[i, ] = foo[3:4]
	}
	
	segments(new.from[, 1], new.from[, 2], new.to[, 1], new.to[, 2])
	
	for(i in 1:ecount(g)) {
		a.coor = arrow.coor(new.from[i, 1], new.from[i, 2],
		                    new.to[i, 1], new.to[i, 2])
		polygon(a.coor, col="black")
	}
	
}



arrow.coor = function(x1, y1, x2, y2, length = 0.02) {
	d = sqrt((x1 - x2)^2 + (y1 - y2)^2)
	st = -(y2 - y1)/d
	ct = (x2 - x1)/d
	A = matrix(c(ct, -st, st, ct), 2)

	M = c(x2, y2)
	p1 = c(-cos(pi/12)*length, sin(pi/12)*length)
	a1 = A %*% p1 + M
	a1 = as.vector(a1)
	p2 = c(-cos(pi/12)*length, -sin(pi/12)*length)
	a2 = A %*% p2 + M
	a2 = as.vector(a2)
	
	return(rbind(a1, a2, c(x2, y2)))
}

reedge = function(x1, y1, x2, y2, r1, r2) {
	d = sqrt((x1 - x2)^2 + (y1 - y2)^2)
	st = -(y2 - y1)/d
	ct = (x2 - x1)/d
	A = matrix(c(ct, -st, st, ct), 2)

	M = c(x1, y1)
	p1 = c(r1, 0)
	a1 = A %*% p1 + M
	a1 = as.vector(a1)
	p2 = c(d - r2, 0)
	a2 = A %*% p2 + M
	a2 = as.vector(a2)
	
	return(c(a1, a2))
}
