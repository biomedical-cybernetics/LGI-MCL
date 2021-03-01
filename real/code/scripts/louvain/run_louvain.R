args = commandArgs(trailingOnly=TRUE)
workdir = args[1]
fileinput = args[2]
fileoutput = args[3]
sparse_flag = strtoi(args[4])

setwd(workdir)
suppressWarnings(suppressMessages(library("igraph")))
suppressWarnings(suppressMessages(library("Matrix")))

# create unweighted matrix
m = as.matrix(read.table(fileinput,header=FALSE,sep=" "))
if (sparse_flag == 1)
{
    x = sparseMatrix(m[,1], m[,2], x=m[,3], symmetric=TRUE)
    g = graph.adjacency(x, mode="undirected", weighted=TRUE, diag=FALSE)
} else
{
    g = graph.adjacency(m, mode="undirected", weighted=TRUE, diag=FALSE)
}

# run louvain and save communities
comm = (multilevel.community(g))$memberships
comm = comm[nrow(comm):1, ]
write.table(comm, file=fileoutput, row.names=FALSE, col.names=FALSE)