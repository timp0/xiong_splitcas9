##Plot results from "methpat" from Bernie Pope (http://bjpop.github.io/methpat/) - to analyze methylation patterns -
##not currently used in paper

library(ggplot2)
library(plyr)
library(reshape2)

library(gtable)
library(grid)

rawdir="/mithril/Data/NGS/Aligned/160720_tina"
plotdir="~/Dropbox/Data/Genetics/Targeted/160803_pattern"

for (k in c("TX20", "TX21", "TX22")) {
    
    for (j in c("top", "bot")) {
        raw=read.delim(file.path(rawdir, paste0(k, ".", j, ".methpat.tsv")))
        
        for (i in unique(raw$X.amplicon.ID.)) {
            
            first=raw[raw$X.amplicon.ID.==i,]
            
            cgpats=data.frame(num=first$X.count., do.call(rbind, (strsplit(as.character(first$X.Methylation.pattern.), split=''))))
            
            
            ##Must have at least 1% of the reads covering this CG (filter indel created CGs)
            thresh=colSums((cgpats[,2:dim(cgpats)[2]]!="-")*cgpats$num)>(.01*sum(cgpats$num))
            
            filt=cgpats[,c(TRUE, thresh)]
            
            ##Filter out reads without all CGs covered
            covpats=rowSums(filt[,2:dim(filt)[2]]!="-")>dim(filt)[2]-2
            
            first=first[covpats,]
            filt=filt[covpats,]
            
            ##Rename columns as CG locs
            locs=strsplit(as.character(first$X.raw.cpg.sites.[1]), split='[(,)]')[[1]]        
            locs=locs[locs>1]
            colnames(filt)[2:dim(filt)[2]]=locs
            
            filt=filt[order(-filt$num),]
            filt$idx=1:dim(filt)[1]
            
            filt.m=melt(filt, id.vars=c("idx", "num"), measure.vars=locs)
            
            pdf(file.path(plotdir, paste0(k, ".", i, ".", j, ".try.pdf")))    
            
            tiley=ggplot(filt.m, aes(y=variable, x=idx, fill=value))+geom_tile(color="black")+theme_bw()+coord_equal()+
                scale_fill_brewer(palette="Set1")
            
            histy=ggplot(filt.m, aes(x=idx, y=num))+geom_bar(stat="identity")+theme_bw()
            
            g1=ggplotGrob(histy)
            g1=gtable_add_cols(g1, unit(0, "mm"))
            g2=ggplotGrob(tiley)
            
            g <- rbind(g1, g2, size="first") # stack the two plots
            g$widths <- unit.pmax(g1$widths, g2$widths) # use the largest widths
            ## center the legend vertically
            g$layout[grepl("guide", g$layout$name),c("t","b")] <- c(1,nrow(g))
            grid.draw(g)
            
            dev.off()
            
            write.table(filt, file=file.path(plotdir, paste0(k, ".", i, ".", j, ".try.csv")), sep=",", row.names=F)
            
        }
    }
    
}




