#'
#'
#'
#'
#'
#'
#'
autoplot.NanoStringRccSet <- function(object, 
    type = c("seq-saturation",
        ), 
    log2scale = TRUE, 
    elt = "exprs", 
    index = 1L, 
    geomParams = list(), 
    tooltipDigits = 4L, 
    heatmapGroup = NULL, 
    blacklist = NULL, 
    tooltipID = NULL, 
    qcCutoffs = list(Housekeeper = c(failingCutoff = 32, 
        passingCutoff = 100), Imaging = c(fovCutoff = 0.75), BindingDensity = c(minimumBD = 0.1, 
        maximumBD = 2.25, maximumBDSprint = 1.8), ERCCLinearity = c(correlationValue = 0.95), 
        ERCCLoD = c(standardDeviations = 2)), scalingFactor = 1L, show_rownames_gene_limit = 60L, 
    show_colnames_gene_limit = 36L, 
    show_rownames_sig_limit = 60L, 
    show_colnames_sig_limit = 36L, 
    subSet = NULL, ...) {

    match.arg()

}

