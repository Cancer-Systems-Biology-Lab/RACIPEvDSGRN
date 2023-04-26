library(tidyverse)

solutionFileProcessMat <- function(net, tol = 1e-3, paramIDs = T, basinThresh = NULL)
{
    ### Adding parameter indices ----
    Net <- net %>% str_remove("\\d") %>% str_remove("_.*")
    # setwd(net)
    nodes <- read.delim(paste0(net, ".prs")) %>% filter(str_detect(Parameter, "Prod")) %>% 
            select(Parameter) %>% unlist %>% str_remove("Prod_of_")
    solnNames <- c("parIndex", "nStatesRAC", "basin", nodes)
    solnFile <- read_csv(paste0(net, "_ss.csv"), col_types = "d")
    colnames(solnFile)[1] <- "parIndex"
    sln <- solnFile %>% group_by(parIndex) %>% summarise(nStatesRAC = n())
    solnFile <- merge(solnFile, sln, by = "parIndex") %>%
        mutate(basin = Count) %>%
        select(all_of(solnNames))
    bMax <- solnFile %>% filter(nStatesRAC == 1) %>% select(basin) %>% unlist %>% unique
    solnFile <- solnFile%>% mutate(basin = basin/bMax)
    pars <- solnFile$parIndex %>% unique
    
    parIds <- data.frame(DSGRNPhase = readLines(paste0(net, "_parameters_FPs.dat")) %>% 
            str_split("\t") %>% 
            sapply(function(x) {
                x %>% unique %>% sort %>% paste0(collapse = "-")
            })) %>%
            mutate(parIndex = 1:nrow(.))
    


    solnFile <- merge(parIds, solnFile, by = "parIndex", all.y = T)
    
    ### discretizing states ----
    nodes <- solnNames[-c(1:3)]
    paramNames <- c("parIndex", "nStatesRAC",
                    read.delim(paste0(net, ".prs")) %>% 
                        select(Parameter) %>% unlist)
    pars <- read_delim(paste0(net, "_parameters.dat"),
                             col_names = paramNames,
                             col_types = "d") 
    parameters <- pars %>% select(parIndex, contains("Trd"))
    links <- colnames(parameters)[-1]   
    OutGoins <- lapply(nodes, function(n){
        links[str_detect(links, paste0(n,"To"))]
    }) %>% set_names(nodes)
    states <- sapply(nodes, function(n){
        sl <- solnFile %>% select(parIndex, all_of(n))
        par <- parameters %>% select(parIndex, all_of(OutGoins[[n]]))
        i <- ncol(par)-1
        d <- merge(sl, par, by = "parIndex", all.x = T) %>% 
            set_names(c("parIndex", "Node", paste0("Link", 1:i)))
        node <- d$Node
        d1 <- d %>% select(contains("Link"))
        rowSums(node > d1)
    }) %>% data.frame %>% set_names(nodes) %>% 
        unite("State", all_of(nodes), sep = "") %>%
        cbind.data.frame(solnFile, .)
    
    # states <- states %>% select(-parNum)
    
    statesUp <- states %>% 
        select(parIndex, State, basin) %>%
        group_by(parIndex, State) %>%
        summarise(basin = sum(basin), .groups = "drop") %>%
        mutate(basinDiscret = basin) %>%
        select(-basin)
    states <- states %>%
        merge(statesUp, by = c("parIndex", "State"), all = T)
    if (!is.null(basinThresh)) {
        states <- states %>% filter(basinDiscret > basinThresh)
    }
    # browser()
    phase <- states %>% 
        group_by(parIndex) %>% 
        summarise(Phase = State %>% unique %>% sort %>% paste0(collapse = "-"), .groups = "drop") 
    phasenStates <- phase %>% mutate(nStates = str_count(Phase, "-") + 1)
    
    states <- merge(phasenStates %>% select(parIndex, nStates, Phase), 
                    states, by = "parIndex", all = T)
    
    # df <- merge(pars %>% select(-nStatesRAC), states, by = c("parIndex"), all = T)
    # # browser()
    # df <- normFunc(paste0(net, ".topo"), df, tol)
    # states <- cbind.data.frame(df, states) %>% mutate(Phase = ifelse(nStatesRAC == 10, "Cycle", Phase))
    write_csv(states, paste0(net, "_solutionFormatterMat.csv"), quote = "none")
    phase <- states %>% select(parIndex, Phase) %>% distinct %>% 
        group_by(Phase) %>% summarise(Count = n(),.groups = "drop") %>%
        mutate(Frequency = Count/sum(Count))
    write_csv(phase, paste0(net, "_phasesMat.csv"), quote = "none")
    # setwd("..")
}

# dirs <- list.dirs(".", full.names = F, recursive = F)
# dirs <- dirs[str_detect(dirs, "GRHL")]
# sapply(dirs, function(d) {
#     nums <- d %>% str_remove(".*_") %>% str_split("To") %>%
#         unlist %>%
#         as.numeric
#     if (all(nums < 60)) {
#         setwd(d)
#         nets <- str_replace(d, "GRHL", c("GRHL", "GRHL1", "GRHL2"))
#         sapply(nets, solutionFileProcessMat)
#         setwd("..")
#     }
# })
