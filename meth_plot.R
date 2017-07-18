##Basic code to analyze and plot methylation data

library(ggplot2)
library(Biostrings)

stdplot <- function(datdir, plotdir, plotreg, samp) {
    ##Make typical plot
    
    cyto.names=c("chr", "pos", "strand", "meth", "unmeth", "context", "trinuc")
    
    files=file.path(datdir, samp, paste0(samp, ".cyto.txt"))
    
    
    mdat=read.delim(files, header=F, col.names=cyto.names)
    mdat$cov=mdat$meth+mdat$unmeth
    mdat$ratio=mdat$meth/mdat$cov

    zmax=max(mdat$ratio[(mdat$pos>plotreg[1])&(mdat$pos<plotreg[2])])
    
    pdf(file.path(plotdir, paste0(samp, ".pdf")), height=8.5, width=11)
    
    print(ggplot(mdat, aes(x=pos, y=ratio, color=context, shape=strand, group=interaction(strand, context)))+theme_bw()+
          geom_point() + geom_line())
    
    print(ggplot(mdat, aes(x=pos, y=ratio, color=context, shape=strand, group=interaction(strand, context)))+theme_bw()+
          geom_point() + geom_line()+xlim(plotreg[1], plotreg[2])+ylim(0, zmax))
    
    
    print(ggplot(mdat, aes(x=pos, y=cov, color=context, shape=strand, group=interaction(strand, context)))+theme_bw()+
          geom_point() + geom_line())
    
    print(ggplot(mdat, aes(x=cov, color=context, shape=strand, group=interaction(strand, context)))+theme_bw()+
          geom_freqpoly()+scale_x_log10())
        
    dev.off()

}

shftplot <- function(datdir, plotdir, plotreg, samp, ref=NA) {
    ##Make plot using "shifted" or pooled plasmid data
    
    if (is.na(ref)) {
        ref=samp
    }
    
    
    cyto.names=c("chr", "pos", "strand", "meth", "unmeth", "context", "trinuc")
    
    files=file.path(datdir, samp, paste0(samp, ".cyto.txt"))
    
    mdat=read.delim(files, header=F, col.names=cyto.names)
    mdat$cov=mdat$meth+mdat$unmeth
    mdat$ratio=mdat$meth/mdat$cov

    seqref=readDNAStringSet(file.path(plotdir, "References", paste0(ref, ".fasta")))

    pat="TGCGCA"
    plotreg=c(1500, 1650)
    cas9reg=IRanges(start=plotreg[1], end=plotreg[2])

    
    cas9loc=unlist(vmatchPattern(pat, seqref))
    cas9loc=subsetByOverlaps(cas9loc, cas9reg)
   

    mdat=mdat[(mdat$pos>plotreg[1])&(mdat$pos<plotreg[2]),]

    ##3rd and 4th position extract

    targs=data.frame()
    
    for (i in 1:length(cas9loc)) {

        plasloc=cas9loc[i]
        
        idx=(names(plasloc)==mdat$chr)&( (mdat$pos==(start(plasloc)+2)) | (mdat$pos==(start(plasloc)+3)))

        subd=mdat[idx,]
        
        subd$nt=i-1
        
        targs=rbind(targs, subd)
        
    }

    
    zmax=max(mdat$ratio[(mdat$pos>plotreg[1])&(mdat$pos<plotreg[2])])
    
    pdf(file.path(plotdir, paste0(samp, ".targonly.pdf")), height=8.5, width=11)

    print(ggplot(targs, aes(x=nt, y=ratio, color=strand))+theme_bw()+
          geom_point() + geom_line())

    print(ggplot(targs, aes(x=nt, y=cov, color=strand))+theme_bw()+
          geom_point() + geom_line())
    
    print(ggplot(targs, aes(x=cov, color=strand))+theme_bw()+
          geom_freqpoly()+scale_x_log10())

    dev.off()

    write.table(targs, row.names=F, file=file.path(plotdir, paste0(samp, ".targonly.tsv")), quote=F, sep="\t")


}


if (FALSE) {

    datdir="/mithril/Data/NGS/Aligned/160325_tina/bismark"
    plotdir="~/Dropbox/Data/Genetics/Targeted/160325_tina"
    
    samp.list=paste0("TX", c(1, 4, 6, 7, 9, 11, 12, 20, 22, 23, 25, 27, 28), "_")

    plotreg=c(850,1850)
    
    for (samp in samp.list) {
        
        stdplot(datdir, plotdir, plotreg, samp)
                
    }
}


if (FALSE) {

    datdir="/mithril/Data/NGS/Aligned/160325_tina/bismark"
    plotdir="~/Dropbox/Data/Genetics/Targeted/160325_tina"

    plotreg=c(1500, 1650)

    samp="TX26_"

    shftplot(datdir, plotdir, plotreg, samp)    
    
}

if (FALSE) {

    datdir="/mithril/Data/NGS/Aligned/160406_tina/bismark"
    plotdir="~/Dropbox/Data/Genetics/Targeted/160407_tina"
    
    samp.list=paste0("TX", c(31, 32), "_")

    plotreg=c(1800, 2800)
    
    for (samp in samp.list) {
        
        stdplot(datdir, plotdir, plotreg, samp)
        
    }

    samp.list=paste0("TX", 33:38, "_")

    plotreg=c(1200, 2200)
    
    for (samp in samp.list) {
        
        stdplot(datdir, plotdir, plotreg, samp)
        
    }

}

if (FALSE) {

    datdir="/mithril/Data/NGS/Aligned/160413_tina/bismark"
    plotdir="~/Dropbox/Data/Genetics/Targeted/160413_tina"


    samp.list=paste0("TX", c(21, 24), "_")
    
    plotreg=c(1200, 2200)
    
    for (samp in samp.list) {
        
        stdplot(datdir, plotdir, plotreg, samp)
        
    }


}

if (FALSE) {

    datdir="/mithril/Data/NGS/Aligned/160622_tina/bismark"
    plotdir="~/Dropbox/Data/Genetics/Targeted/160622_tina"

    samp.list=paste0("TX", 42:47)
    
    plotreg=c(1200, 2200)
    
    for (samp in samp.list) {
        
        stdplot(datdir, plotdir, plotreg, samp)
        
    }

}

if (FALSE) {
    
    datdir="/mithril/Data/NGS/Aligned/160622_tina/bismark"
    plotdir="~/Dropbox/Data/Genetics/Targeted/160622_tina"
    
    plotreg=c(1500, 1650)

    samp="TX41"

    shftplot(datdir, plotdir, plotreg, samp)    

}


if (FALSE) {

    #TXM50- Here we’re testing a new guideRNA. Sample contains one template. Could we get the plot for methylation versus position?

    datdir="/mithril/Data/NGS/Aligned/160728_tina/bismark"
    plotdir="~/Dropbox/Data/Genetics/Targeted/160902_tina"

    samp="tx50"
    
    plotreg=c(1200, 2200)
    
    stdplot(datdir, plotdir, plotreg, samp)
        
    

}

