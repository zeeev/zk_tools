#usage:  nohup R --vanilla < plotLD.R --args my binary distance matrix &
#columns should be indviduals and rows should be loci.  Data should be presence
#absence...

cmd_args <- commandArgs(trailingOnly = TRUE)

plotTREE<-function(x){


	options(error=dump.frames)
	require(ape)                

	dat<-t(read.csv(x, sep="\t", header=TRUE))
	my.dist<-dist(dat)
	my.tree<-nj(my.dist)
	my.boot<-boot.phylo(phy=my.tree, dat, FUN=function(xx){nj(dist(xx))}, B=100, quite=FALSE, check.labels=FALSE)
	write.csv(file="boot.vales.txt", my.boot)		
	write.tree(file="tree.out", my.tree)
}


plotTREE(cmd_args)

