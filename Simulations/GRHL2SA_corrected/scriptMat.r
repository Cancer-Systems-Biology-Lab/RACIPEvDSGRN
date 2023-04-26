source("solutionFormatMatlab.r")
source("compareStabilityMat.R")
library(funcsKishore)

dirs <- list.dirs(".", full.names = F, recursive = F)
dirs <- dirs[str_detect(dirs, "GRHL")]
sapply(dirs, function(d) {
    nums <- d %>% str_remove(".*_") %>% str_split("To") %>%
        unlist %>%
        as.numeric
    if (all(nums < 60)) {
        print(d)
        setwd(d)
        nets <- str_replace(d, "GRHL", c("GRHL", "GRHL1", "GRHL2"))
        sapply(nets, solutionFileProcessMat)
        sapply(nets, compareMat)
        setwd("..")
        # jsdCompiled
    }
})

params <- c("Phase", "monoMulti", "nState")
nStates <- 1:4
nList1 <- list(c(1), 1:3, 2:3, 2:4, 3:5, 4:6, 5:7, 6:8, 7:9, 8:10, 9:11, 10:12)
nList2 <- lapply(seq(2, 50, 5), function(i) {
    i:(i+9)
})
nList3 <- lapply(seq(1, 41, 10), function(i) {
    i:(i+9)
})


nLists <- list(nList1, nList2, nList3, c(nList1, nList2, nList3))
names(nLists) <- c("SmallMat", "MediumMat", "LargeMat", "AllMat")

zeroNa <- function(d) {
    d[is.na(df)] <- 0
    d
}

sapply(names(nLists), function(nListName) {
    # browser()
    nList <- nLists[[nListName]]
    nets <- sapply(nList, function(x) {
        paste0("GRHL_", x[1], "To", x[length(x)])
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
                key <- str_remove(net, "GRHL_")
                keys <- paste0(c("GRHL","GRHL1", "GRHL2"), "_",key)
                jsdFils <- paste0(keys, "_JSDdataMat.csv")
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
    setwd("..")
    write_csv(df, paste0("GRHL_", nListName, "_allData.csv"))    
})
