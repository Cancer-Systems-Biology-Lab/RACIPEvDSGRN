library(tidyverse)
library(funcsKishore)
params <- c("Phase", "monoMulti", "nState")
nStates <- 1:4
nList1 <- list(c(1), 1:3, 2:3, 2:4, 3:5, 4:6, 5:7, 6:8, 7:9, 8:10, 9:11, 10:12)
nList2 <- lapply(seq(2, 55, 5), function(i) {
    i:(i+9)
})
nList3 <- lapply(seq(1, 41, 10), function(i) {
    i:(i+9)
})

nLists <- list(nList1, nList3, c(nList1, nList3))
names(nLists) <- c("Small", "Large", "All")

zeroNa <- function(d) {
    d[is.na(df)] <- 0
    d
}

sapply(names(nLists), function(nListName) {
    # browser()
    nList <- nLists[[nListName]]
    nets <- sapply(nList, function(x) {
        paste0("TSq_", x[1], "To", x[length(x)])
    })
    xLabs <- sapply(nList, function(x) {
        if (length(x) == 1) {
            as.character(x)
        }
        else {
            paste0(x[1], "-", x[length(x)])
        }
    })
    df  <- lapply(1:length(nList), function(i) {
                nL <- nList[[i]]
                net <- nets[[i]]
                xLab <- xLabs[[i]]
                setwd(net)
                key <- str_remove(net, "TSq_")
                keys <- paste0(c("TSq","TSq1", "TSq2"), "_",key)
                jsdFils <- paste0(keys, "_JSDdata.csv")
                jsdCompiled <- lapply(jsdFils, function(j) {
                    read_csv(j) %>% select(-RACIPE, -DSGRN, -Phase) %>% unique
                }) %>% reduce(merge, by = c("Parameter","nState"), all = T)
                jsdCompiled <- jsdCompiled %>% 
                    # zeroNa %>%
                    mutate(Avg = jsdCompiled %>% select(contains("JSD")) %>% apply(1, mean),
                        Std = jsdCompiled %>% select(contains("JSD")) %>% apply(1, sd),
                        HillRange = xLab) %>%
                    select(-contains("JSD"))
                setwd("..")
                jsdCompiled
            }) %>%
            reduce(rbind.data.frame) %>%
            mutate(HillRange = factor(HillRange, levels = xLabs), 
                nState = as.character(nState))
    dfnState <- lapply(1:length(nList), function(i) {
        nL <- nList[[i]]
        net <- nets[[i]]
        xLab <- xLabs[[i]]
        setwd(net)
        key <- str_remove(net, "TSq_")
        keys <- paste0(c("TSq","TSq1", "TSq2"), "_",key)
        FPs <- paste0(keys, "_parameters_FPs.dat")
        FPDat <- lapply(FPs, function(f) {
            readLines(f) %>% str_split("\t") %>% sapply(length) %>%
            table %>% data.frame %>% set_names(c("nStates", "Freq")) %>%
                mutate(Freq = Freq/sum(Freq))
        }) %>% reduce(merge, by = "nStates", all = T) %>% zeroNa
        setwd("..")
        FPDat <- FPDat %>%
            mutate(Avg = FPDat %>% select(contains("Freq")) %>% apply(1, mean),
                   Std = FPDat %>% select(contains("Freq")) %>% apply(1, sd),
                   HillRange = xLab)
    }) %>% reduce(rbind.data.frame) %>% 
        mutate(HillRange = factor(HillRange, levels = xLabs), 
               nStates = as.character(nStates))
    dfPhaseDS <- lapply(1:length(nList), function(i) {
        nL <- nList[[i]]
        net <- nets[[i]]
        xLab <- xLabs[[i]]
        setwd(net)
        key <- str_remove(net, "TSq_")
        keys <- paste0(c("TSq","TSq1", "TSq2"), "_",key)
        FPs <- paste0(keys, "_parameters_FPs.dat")
        FPDat <- lapply(FPs, function(f) {
            readLines(f) %>% str_split("\t") %>% sapply(function(x) {
                x %>% unique %>% sort %>% paste0(collapse = "-")
            }) %>%
                table %>% data.frame %>% set_names(c("Phase", "Freq")) %>%
                mutate(Freq = Freq/sum(Freq))
        }) %>% reduce(merge, by = "Phase", all = T) %>% zeroNa
        setwd("..")
        FPDat <- FPDat %>%
            mutate(Avg = FPDat %>% select(contains("Freq")) %>% apply(1, mean),
                   Std = FPDat %>% select(contains("Freq")) %>% apply(1, sd),
                   HillRange = xLab)
    }) %>% reduce(rbind.data.frame) %>% 
        mutate(HillRange = factor(HillRange, levels = xLabs), 
               nStates = str_count(Phase, "-") + 1)
    write_csv(dfPhaseDS, paste0(nListName, "_PhaseDat.csv"))
    dir.create(nListName)
    setwd(nListName)
    sapply(params, function(p) {
        sapply(nStates, function(nSt) {
            dfFil <- df %>% filter(Parameter == p, nState == nSt)
            ggplot(dfFil, aes(x = HillRange, y = Avg)) +
                geom_bar(stat = "identity") +
                geom_errorbar(aes(ymin = Avg-Std, ymax = Avg+Std)) +
                theme_Publication() + labs(x = "Hill Range", y = "JSD") +
                theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
            ggsave(paste0(nListName, "_", p, "_", nSt, ".png"), width = 5.5, height = 5)
        })
        dfFil <- df %>% filter(Parameter == p)
        ggplot(dfFil %>% filter(nState %in% nStates), aes(x = HillRange, y = Avg, fill = nState)) +
                geom_bar(stat = "identity", position = position_dodge()) +
                geom_errorbar(aes(ymin = Avg-Std, ymax = Avg+Std), position = position_dodge(width = 0.9)) +
                theme_Publication() + labs(x = "Hill Range", y = "JSD") +
                theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
        ggsave(paste0(nListName, "_", p, "_allN", ".png"), width = 5.5, height = 5)
    })
    ggplot(dfnState, aes(x = HillRange, y = Avg, fill = nStates)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        geom_errorbar(aes(ymin = Avg-Std, ymax = Avg + Std), position = position_dodge(width = 0.9)) +
        theme_Publication() + labs(x = "Hill Range", y = "Frequency") +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
    ggsave(paste0(nListName, "_nStateFreq", ".png"), width = 5.5, height = 5)
    setwd("..")
    write_csv(df, paste0("TSq_", nListName, "_allData.csv"))    
})
