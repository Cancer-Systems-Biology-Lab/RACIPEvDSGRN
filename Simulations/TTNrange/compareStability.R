library(philentropy)
# library(tidyverse)
library(funcsKishore)

# nets <- list.files(".", ".topo") %>% str_remove(".topo")
# solutionFiles <- lapply(nets, function(net) {
#     paste0(net, "_solutionFormatter.csv") %>% read_csv
# })
# names(solutionFiles)  <- nets

# parameterFiles <- lapply(nets, function(net) {
#     paramNames <- paste0(net, ".prs") %>% read.delim %>% 
#         select(Parameter) %>% unlist %>% c("parIndex", "nStatesRAC", .)
#     paste0(net, "_parameters.dat") %>% read_delim(delim = "\t", col_types = "d", 
#         col_names = paramNames)
# })
# names(parameterFiles) <- nets
compare <- function(net, basinThresh = NULL) {
    params <- c("Phase", "monoMulti", "nState")
    solnFile <- paste0(net, "_solutionFormatter.csv") %>%
        read_csv %>%
        mutate(nStatesRAC = str_count(Phase, "-") + 1,
        nStatesDSGRN = str_count(DSGRNPhase, "-") + 1,
        basin = basin/max(basin)) %>%
        mutate(stRACIPE = ifelse(nStatesRAC == 1, "Monostable", "Multistable"),
            stDSGRN = ifelse(nStatesDSGRN == 1, "Monostable", "Multistable"))
    if (!is.null(basinThresh) && is.numeric(basinThresh)) {
        solnFile <- solnFile %>% 
            filter(basin > basinThresh) %>% select(-Phase)
        phases <- solnFile %>% group_by(parIndex) %>% summarise(Phase = State %>% unique %>% sort %>% paste0(collapse = "-"))
        solnFile <- merge(solnFile, phases, by = "parIndex")

    }
    nStateList <- c(unique(solnFile$nStatesDSGRN), NULL)
    fDf <- lapply(nStateList, function(nState) {
        if (!is.null(nState)) {
            solnFile <- solnFile %>% filter(nStatesDSGRN == nState)
        }
        else {
            nState = "All"
        }
        lapply(params, function(param) {
            if (param == "Phase") {
                racCol <- "Phase"
                DSGCol <- "DSGRNPhase"
            }
            else if (param == "monoMulti") {
                racCol <- "stRACIPE"
                DSGCol <- "stDSGRN"
            }
            else {
                racCol <- "nStatesRAC"
                DSGCol <- "nStatesDSGRN"
            }
            colnames(solnFile)[colnames(solnFile) == racCol]  <- "racCol"
            colnames(solnFile)[colnames(solnFile) == DSGCol]  <- "DSGCol"
            # browser()
            rac <- solnFile %>% select(parIndex, racCol) %>% 
                unique %>% group_by(racCol) %>% summarise(RACIPE = n()) %>%
                mutate(RACIPE = RACIPE/sum(RACIPE)) %>%
                set_names(c("Phase", "RACIPE"))
            dsgrn <- solnFile %>% select(parIndex, DSGCol) %>%
                unique %>% group_by(DSGCol) %>% summarise(DSGRN = n()) %>%
                mutate(DSGRN = DSGRN/sum(DSGRN)) %>% 
                set_names(c("Phase", "DSGRN"))
            phaseDf <- merge(rac, dsgrn, by = "Phase", all = T)
            phaseDf[is.na(phaseDf)] <- 0
            phaseDf %>% 
                mutate(JSD = JSD(phaseDf %>% select(-Phase) %>% t),
                    Parameter = param, nState = nState)
        }) %>% reduce(rbind.data.frame) 
        # phaseDf <- pDf
    }) %>% reduce(rbind.data.frame)
    phaseDf <- fDf
    
    write_csv(phaseDf, paste0(net, "_JSDdata.csv"))
}
