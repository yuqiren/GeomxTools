DEFAULTS <- list(minSaturation=0.7, minReads=10000, minProbeRatio=0.1, 
    minimumCount=10, localOutlierAlpha=0.01, globalOutlierRatio=0.2, 
    loqMultiplier=2.0, loqCutoff=0.7, highCountCutoff=10000)

#' Add QC flags to feature or protocol data
#' 
#' @param object name of the object class to perform QC on
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
#' }
#' @param qcCutoffs list of cutoffs and thresholds to use for QC
#' 
#' @return the object that QC was performed on
#' 
#' @examples
#' 
setMethod("setQCFlags",
    signature(object="NanoStringGeoMxSet"),
    function(object, qcCutoffs=DEFAULTS, ...) {
        qcCutoffs <- checkCutoffs(qcCutoffs)
        object <- setAOIFlags(object=object, qcCutoffs=qcCutoffs)
        if (featureType(object) == "Probe") {
            object <- setProbeFlags(object=object, qcCutoffs=qcCutoffs)
        } else if (featureType(object) == "Target") {
            object <- setTargetFlags(object=object, qcCutoffs=qcCutoffs)
        } else {
            valid(Object(x))
        }
        return(object)
})

setAOIFlags <- function(object, qcCutoffs=DEFAULTS) {
    object <- setSaturationFlags(object=object, 
        cutoff=qcCutoffs[["minSaturation"]])
    object <- setLowReadFlags(object=object, 
        cutoff=qcCutoffs[["minReads"]])
    return(object)
}

setProbeFlags <- function(object, qcCutoffs=DEFAULTS) {
    object <- setProbeRatioFlags(object=object, 
        cutoff=qcCutoffs[["minProbeRatio"]])
    object <- setProbeCountFlags(object=object, 
        cutoff=qcCutoffs[["minimumCount"]])
    object <- setLocalFlags(object=object, 
        cutoff=qcCutoffs[["localOutlierAlpha"]])
    object <- setGlobalFlags(object=object, 
        cutoff=qcCutoffs[["globalOutlierRatio"]])
    return(object)
}

setTargetFlags <- function(object, qcCutoffs=DEFAULTS) {
    object <- 
        setLOQFlags(object=object, multiplier=DEFAULTS[["loqMultiplier"]], 
            cutoff=qcCutoffs[["loqCutoff"]])
    #object <- 
    #    setHighCountFlags(object=object, cutoff=qcCutoffs[["highCountCutoff"]])
    return(object)
}

setSaturationFlags <- function(object, cutoff=DEFAULTS[["minSaturation"]]) {
    percentUnique <- 
        sData(object)["DeduplicatedReads"] / sData(object)["Aligned"]
    percentUnique <- percentUnique > cutoff
    colnames(percentUnique) <- "Saturation"
    object<- appendSampleFlags(object, percentUnique)
    return(object)
}

setLowReadFlags <- function(object, cutoff=DEFAULTS[["minReads"]]) {
    lowReads <- sData(object)["Raw"] < cutoff
    colnames(lowReads) <- "LowReads"
    object <- appendSampleFlags(object, lowReads)
    return(object)
}

setProbeRatioFlags <- 
    function(object=object, cutoff=DEFAULTS[["minProbeRatio"]]) {
        rawTargetCounts <- collapseCounts(object)
        rawTargetCounts[["Mean"]] <- 
            apply(rawTargetCounts[, sampleNames(object)], 
                MARGIN=1, FUN=ngeoMean)
        rownames(rawTargetCounts) <- rawTargetCounts[["TargetName"]]
        targetMeans <- rawTargetCounts[fData(object)[["TargetName"]], "Mean"]
        probeMeans <- apply(assayDataElement(object, elt="exprs"), 
            MARGIN=1, FUN=ngeoMean)
        probeRatioFlags <- (probeMeans / targetMeans) < cutoff
        probeRatioFlags <- data.frame("LowProbeRatio"=probeRatioFlags)
        object <- appendFeatureFlags(object, probeRatioFlags)
        return(object)
    }
 
setProbeCountFlags <- 
    function(object=object, cutoff=DEFAULTS[["minimumCount"]]) {
        probeCountFlags <- apply(assayDataElement(object, elt="exprs"), 
            MARGIN=1, FUN=function(x, minCount){
                all(x < minCount)
            }, minCount=cutoff)
        probeCountFlags <- data.frame("LowProbeCount"=probeCountFlags)
        object <- appendFeatureFlags(object, probeCountFlags)
        return(object)
    }

setLocalFlags <- 
    function(object=object, cutoff=DEFAULTS[["localOutlierAlpha"]]) {
        #if ("LowProbeCount" %in% names(fData(object)[["QCFlags"]])) {
          #  subObj <- object[!fData(object)[["QCFlags"]][["LowProbeRatio"]], ]
            probeCounts <- 
                setDT(cbind(fData(object)[, c("RTS_ID", "TargetName", "Module")], 
                    assayDataElement(object, elt="exprs")))
            probeCounts <- melt(probeCounts, 
                id.vars=c("RTS_ID", "TargetName", "Module"), 
                variable.name="Sample_ID", 
                value.name="Count", variable.factor=FALSE)
            probeCounts[, Count:=logt(Count)]
            probeCounts[, "LowLocalOutlier"] <- FALSE
            probeCounts[, "HighLocalOutlier"] <- FALSE
            probeCounts <- probeCounts[, suppressWarnings(grubbsFlag(.SD, alpha=cutoff)), 
                by=.(TargetName, Module, Sample_ID)]
            lowFlags <- as.data.frame(dcast(probeCounts, RTS_ID ~ Sample_ID, value.var="LowLocalOutlier"), stringsAsFactor=FALSE)
            highFlags <- as.data.frame(dcast(probeCounts, RTS_ID ~ Sample_ID, value.var="HighLocalOutlier"), stringsAsFactor=FALSE)
            rownames(lowFlags) <- lowFlags[["RTS_ID"]]
            rownames(highFlags) <- highFlags[["RTS_ID"]]
            lowFlags <- lowFlags[, colnames(lowFlags) != "RTS_ID"]
            highFlags <- highFlags[, colnames(highFlags) != "RTS_ID"]
            outlierFlags <- data.frame(LowLocalOutlier=lowFlags, HighLocalOutlier=highFlags)
            object <- appendFeatureFlags(object, outlierFlags)
        #} #else {
        #    stop(paste("It is recommended to flag targets with overall", 
        #        "low counts (i.e. background-level expression)", 
        #        "prior to checking for outliers.\n", "Rerun outlier testing",
        #        "after running setProbeCountFlags."))
        #}
        return(object)
    }

setGlobalFlags <- 
    function(object=object, cutoff=DEFAULTS[["globalOutlierRatio"]]) {
        lowFlagRatio <- 
            apply(fData(object)[["QCFlags"]][, grepl("LowLocalOutlier", 
                colnames(fData(object)[["QCFlags"]]))], 1, mean )
        highFlagRatio <- 
            apply(fData(object)[["QCFlags"]][, grepl("HighLocalOutlier", 
                colnames(fData(object)[["QCFlags"]]))], 1, mean )
        outlierRatio <- lowFlagRatio > cutoff | highFlagRatio > cutoff
        globalFlags <- data.frame("GlobalOutlier"=outlierRatio)
        object <- appendFeatureFlags(object, globalFlags)
        return(object)
    }

setLOQFlags <- 
    function(object=object, multiplier=DEFAULTS[["loqMultiplier"]], cutoff=DEFAULTS[["loqCutoff"]]) {
        if (featureType(object) == "Target") {
            LOQs <- 
                pData(object)[, "NegGeoMean", drop=FALSE] * 
                pData(object)[, "NegGeoSD", drop=FALSE] ^ multiplier
            pData(object)[, "LOQ"] <- LOQs
            targetFlags <- t(esApply(object, MARGIN=1, FUN=function(x) {
                return(x < LOQs[sampleNames(object), ])}))
            LOQFlags <- as.data.frame(apply(targetFlags, 1, 
                function(x) {return(mean(x) > cutoff)}))
            colnames(LOQFlags) <- "LOQFlags"
            LOQFlags <- LOQFlags[featureNames(object), , drop=FALSE]
            object <- appendFeatureFlags(object, LOQFlags)
        } else {
            warning(paste("Incorrect feature type.",
                "Feature type should be Target."))
        }
        return(object)
    }

setHighCountFlags <- 
    function(object, cutoff=qcCutoffs[["highCountCutoff"]]) {
}

checkCutoffs <- function(qcCutoffs) {
    if (!all(names(DEFAULTS) %in% names(qcCutoffs))) {
        qcCutoffs <- append(qcCutoffs, 
            DEFAULTS[!(names(DEFAULTS) %in% names(qcCutoffs))])
    }
    return(qcCutoffs)
}

appendSampleFlags <- function(object, currFlags) {
    if("QCFlags" %in% varLabels(protocolData(object))) {
        protocolData(object)[["QCFlags"]] <- 
            cbind(protocolData(object)[["QCFlags"]], currFlags) 
    } else {
        protocolData(object)[["QCFlags"]] <- currFlags
    }
    return(object)
}

appendFeatureFlags <- function(object, currFlags) {
    if("QCFlags" %in% varLabels(featureData(object))) {
        featureData(object)[["QCFlags"]] <- 
            cbind(featureData(object)[["QCFlags"]], currFlags) 
    } else {
        featureData(object)[["QCFlags"]] <- currFlags
    }
    return(object)
}

grubbsFlag <- function(countDT, alpha=0.01) {
    if (dim(countDT)[1] < 3 | all(countDT[, Count] == countDT[1][, Count])) {
        return(countDT)
    }
    grubbsResult <- outliers::grubbs.test(countDT[, Count], two.sided=TRUE)
    if (grubbsResult$p.value < alpha) {
        if (grepl("lowest", tolower(grubbsResult$alternative))) {
            countDT[Count == min(Count), "LowLocalOutlier"] <- TRUE
        } else {
            countDT[Count == max(Count), "HighLocalOutlier"] <- TRUE
        }
    } else {
        return(countDT)
    }
    return(countDT)
}
