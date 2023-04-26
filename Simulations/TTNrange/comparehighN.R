library(philentropy)
library(tidyverse)
library(funcsKishore)

setwd("D:\\Github\\Projects\\Ongoing\\DSGRN\\cleanVersion/ParameterSets/TT")
FPList <- lapply(c("TT", "TT1", "TT2"), function(x) {
    readLines(paste0(x, "_parameters_FPs.dat"))
})
names(FPList) <- c("TT", "TT1", "TT2")
setwd("D:\\Github\\Projects\\Ongoing\\DSGRN\\cleanVersion/SmallSample/TT")

compare <- function(net, n1, n2, n3, n4, n5, n6, 
    nState = NULL, param = c("Phase", "monoMulti", "nState")) {
        # browser()
    param  <- match.arg(param)
    Net <- net %>% str_extract("TT.?_") %>% str_remove("_")
    Fps <- FPList[[Net]]
    paramNames <- paste0(net, ".prs") %>% read.delim %>% 
        select(Parameter) %>% unlist %>% c("parIndex", "nStatesRAC", .)
    parIDs <- paste0(net, "_parameters.dat") %>%
        read_delim(delim = "\t", col_types =  "d", col_names = paramNames) %>%
        filter(Num_of_AToB %in% n1, Num_of_BToA %in% n2,
                Num_of_AToC %in% n3, Num_of_CToA %in% n4,
                Num_of_BToC %in% n5, Num_of_CToB %in% n6) %>% 
        select(parIndex) %>% unlist
        # browser()
    solnFile <- paste0(net, "_solutionFormatter.csv") %>%
        read_csv(col_types = cols()) %>%
        filter(parIndex %in% parIDs) %>%
        mutate(DSGRNPhase = Fps[parIndex] %>% str_replace_all("\t", "-")) %>%
        mutate(nStatesRAC = str_count(Phase, "-") + 1,
        nStatesDSGRN = str_count(DSGRNPhase, "-") + 1) %>%
        mutate(stRACIPE = ifelse(nStatesRAC == 1, "Monostable", "Multistable"),
            stDSGRN = ifelse(nStatesDSGRN == 1, "Monostable", "Multistable"))
    if (!is.null(nState)) {
        solnFile  <- solnFile %>% filter(nStatesDSGRN %in% nState)
    }
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
    write_csv(phaseDf, paste0(net, "/",param, "_", nState, paste0("_", n1[1], "-", n1[length(n1)]),".csv"))
    JSD(phaseDf %>% select(-Phase) %>% t)
}

sapply(c(1:3, NULL), function(nState) {
    nList1 <- lapply(1:8, function(i) {
        c(i, i+1, i+2)
    })

    nList2 <- lapply(seq(11,40, 5), function(i) {
        (i):(i+9)
    })
    sapply(c("Phase", "nState", "monoMulti"), function(param) {
        df <- lapply(c("TT", "TT1", "TT2"), function(net) {
            # browser()
            Net <- paste0(net, "_1To10")
            dir.create(Net, showWarnings = F)
            l <- sapply(nList1, function(n) {
                compare(Net, n,n,n,n,n,n, nState, param)
            })
            data.frame(nRange = sapply(nList1, function(k) {
                paste0(k[1], "-", k[length(k)])
            }),
                JSD = l) %>% set_names(c("nRange", net))
        }) %>% reduce(merge, by = "nRange", all = T)

        df <- df %>% 
            mutate(Avg = df %>% select(-nRange) %>% apply(1, mean),
                Std = df %>% select(-nRange) %>% apply(1, sd))

        ggplot(df, aes(x = nRange, y = Avg)) +
            geom_bar(stat = "identity") +
            geom_errorbar(aes(ymin = Avg-Std, ymax = Avg+Std)) +
            theme_Publication() + 
            theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
            labs(x = "Range of Hill Coefficients", y = "JSD")

        ggsave(paste0("TT1To10_JSD_", param, nState,".png"), width = 5.5, height = 5.5)
    })
    sapply(c("Phase", "nState", "monoMulti"), function(param) {
        df <- lapply(c("TT", "TT1", "TT2"), function(net) {
            Net <- paste0(net, "_10To50")
            dir.create(Net, showWarnings = F)
            l <- sapply(nList2, function(n) {
                compare(Net, n,n,n,n,n,n, nState, param)
            })
            data.frame(nRange = sapply(nList2, function(k) {
                paste0(k[1], "-", k[length(k)])
            }),
                JSD = l) %>% set_names(c("nRange", net))
        }) %>% reduce(merge, by = "nRange", all = T)

        df <- df %>% 
            mutate(Avg = df %>% select(-nRange) %>% apply(1, mean),
                Std = df %>% select(-nRange) %>% apply(1, sd))

        ggplot(df, aes(x = nRange, y = Avg)) +
            geom_bar(stat = "identity") +
            geom_errorbar(aes(ymin = Avg-Std, ymax = Avg+Std)) +
            theme_Publication() + 
            theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
            labs(x = "Range of Hill Coefficients", y = "JSD")

        ggsave(paste0("TT10To50_JSD_", param, nState,".png"), width = 5.5, height = 5.5)
    })
    
})
