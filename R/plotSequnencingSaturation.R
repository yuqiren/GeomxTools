plotSeqSaturation <-  function(object, 
                               xAxisVar="area", 
                               colorVar="segment", 
                               cutoff=0.7) {
    if (is.numeric(sData(object)[[xAxisVar]])) {
        satPlot <- 
            ggplot(sData(object), mapping=aes_string(x=xAxisVar, 
                y="SequencingSaturation", col=as.factor(colorVar))) +
            geom_point(size=3, alpha=0.4) +
            geom_hline(yintercept=cutoff, size=1, color="hotpink", 
                lty="21", alpha=0.5) +
            ggtitle("Sequencing Saturation", 
                subtitle = paste("By", x=xAxisVar)) +
            scale_x_continuous(labels=function(x) {
                    formatC(x, format="d", big.mark=",")
                },
                limits=c(min(sData(object)[[xAxisVar]]) * 0.8, 
                    max(sData(object)[[xAxisVar]]) * 1.2)) +
            scale_y_continuous(labels=function(x) {
                    paste0(x, "%")
                },
                limits=c(0, 100)) +
            theme_minimal(base_size=16) +
            theme(panel.grid=element_blank(),
                axis.line.x=element_line(),
                axis.line.y=element_line()) +
            labs(x=xAxisVar, y="Sequencing saturation", 
                color="Grouping\nVariable")
        if ("SequencingSaturation" %in% colnames(sData(object)[["QCFlags"]])) {
            satPlot <- satPlot +
                geom_text_repel(data=sData(object[, 
                    !sData(object)[["QCFlags"]][, "SequencingSaturation"]]),
                aes(label=dimLabels(object)[1L]), color="darkgray", 
                box.padding=.25, point.padding=.5, min.segment.length=.5, 
                fontface="bold", direction="y")
        }
        return(satPlot)
    } else {
        satPlot <- ggplot(sData(object), mapping=aes_string(x=xAxisVar,
                y="SequencingSaturation", fill=as.factor(colorVar))) +
            geom_violin(show.legend=FALSE) +
            geom_boxplot(width=0.15, outlier.colour=NA, show.legend=FALSE) +
            geom_hline(yintercept=cutoff, size=1, color="hotpink", 
                lty="21", alpha=0.5) +
            ylim(0,100) +
            ggtitle(paste0("Sequencing Saturation by ", xAxisVar)) +
            theme_minimal(base_size=12) +
            theme(panel.grid=element_blank(), axis.line.x=element_line(),
                axis.line.y=element_line()) +
            labs(x=xAxisVar, y="Sequencing saturation", fill="")
        return(satPlot)
    }
}