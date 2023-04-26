library(funcsKishore)
exact <- T

nList1 <- list(c(1), 1:3, 2:3, 2:4, 3:5, 4:6, 5:7, 6:8, 7:9, 8:10, 9:11, 10:12)
nList2 <- lapply(seq(2, 55, 5), function(i) {
    i:(i+9)
})
nList3 <- lapply(seq(1, 41, 10), function(i) {
    i:(i+9)
})

nLists <- list(nList1, nList2, nList3, c(nList1, nList2, nList3))
names(nLists) <- c("Small", "Medium", "Large", "All")

compareStates <- function(s1, s2) {
    # browser()
    s1 <- s1 %>% str_split("") %>% unlist %>% as.integer
    s2 <- s2 %>% str_split("") %>% unlist %>% as.integer
    abs(s1-s2) %>% sum
}

comparePhase <- function(pL, exact = F) {
    p1 <- pL[1]
    p2 <- pL[2]
    p1 <- str_split(p1, "-") %>% unlist
    p2 <- str_split(p2, "-") %>% unlist
    dist <- sapply(p1, function(p) {
        if (exact) {
            return(p %in% p2)
        }
        l <- sapply(p2, function(q) {compareStates(p,q)}) %>% min
        return(l)
    }) %>% mean
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
    
    dd  <- lapply(1:length(nList), function(i) {
        nL <- nList[[i]]
        net <- nets[[i]]
        xLab <- xLabs[[i]]
        setwd(net)
        key <- str_remove(net, "GRHL_")
        keys <- paste0(c("GRHL","GRHL1", "GRHL2"), "_",key)
        nKey <- "phaseDist"
        if (exact) nKey <- "phaseDistExact"
        racSol <- paste0(keys, "_solutionFormatter.csv") %>% lapply(function(x) {
            df <- read_csv(x) %>% mutate(
                nStateRac = str_count(Phase, "-") +1,
                nStateDSGRN = as.character(str_count(DSGRNPhase, "-")+ 1)
            ) %>% 
                select(parIndex, Phase, DSGRNPhase, nStateRac, nStateDSGRN) %>% unique
            df$Compare <- df %>% select(Phase, DSGRNPhase) %>% t %>% data.frame %>% sapply(function(o) {comparePhase(o, exact)})
            ggplot(df, aes(x = nStateDSGRN, y = Compare)) + geom_violin() +
                theme_Publication() + labs(x = "# states DSGRN", y = "Mean Difference in Phase")
            ggsave(str_replace(x, "_solutionFormatter.csv", paste0("_", nKey, ".png")), width = 5.5, height = 5)
            ggplot(df, aes(x = nStateDSGRN, y = Compare)) + geom_boxplot() +
                theme_Publication() + labs(x = "# states DSGRN", y = "Mean Difference in Phase")
            ggsave(str_replace(x, "_solutionFormatter.csv", paste0("_", nKey, "Box.png")), width = 5.5, height = 5)
            write_csv(df, str_replace(x, "_solutionFormatter.csv", "_phaseDist.csv"))
        })
        setwd("..")
        
    })
})
