library(tidyverse)
# library(funcsKishore)
# options(lazy = F)

wd <- "D:\\Github\\Projects\\Ongoing\\DSGRN\\cleanVersion\\SmallSample"

hill <- function(A, A0, nA, lA)
{
    hminus <- A0^nA/(A^nA + A0^nA)
    hminus + lA*(1-hminus)
}

normFunc <- function(topoFile, df, tol = 1e-3) {
    topoDf <- read.delim(topoFile %>% str_remove("_switching"), sep = "")
    nodes <- unique(c(topoDf$Source, topoDf$Target)) %>% sort
    eqns <- sapply(nodes, function(x) {
        sources <- topoDf %>% filter(Target == x)
        prod <- df[[paste0("Prod_of_", x)]]
        if (length(sources) != 0) {
            apply(sources, 1, function(y) {
                # browser()
                thr <- paste0("Trd_of_", y[1], "To", y[2])
                n <- paste0("Num_of_", y[1], "To", y[2])
                l <- paste0(ifelse(y[3] == "1", "Act_of_", "Inh_of_"), y[1], "To", y[2])
                prod <<- prod*hill(2^(df[[y[1]]]), df[[thr]], df[[n]], df[[l]])
                if (y[3] == "1")
                    prod <<- prod/df[[l]]
            })
        }
        deg <- df[[paste0("Deg_of_", x)]]*(2^(df[[x]]))
        prod - deg
    }) %>% data.frame %>% set_names(paste0("Deriv_",nodes)) %>%
        mutate(Norm = apply(., 1, norm, type = "2")) %>%
        mutate(Steady = ifelse(Norm < tol, "Yes", "No"))
    return(eqns)
}

DSGRNParams <- function(net) {
    browser()
    parFile <- paste0(net, "_parameters.dat")
    if (!file.exists(parFile))
        return()
    prsFile <- paste0(net, ".prs") %>% read.delim
    paramNames <- prsFile %>% select(Parameter) %>% unlist %>% c("parIndex", "nStates", .)
    parDat <- parFile %>% 
        read_delim(delim = "\t", col_names = paramNames, col_types = "d") %>%
        mutate(across(.cols = contains("Act"), .fns = function(x){1/x}))
    colnames(parDat) <- str_replace(colnames(parDat), "Act", "Inh")
    topoFile <- read.delim(paste0(net, ".topo"), sep = "")
    edges <- topoFile %>% 
        select(-Type) %>% 
        mutate(Edge = paste0(Source, "_", Target)) %>%
        select(Edge) %>% unlist
    d <- sapply(edges, function(x){
        y <- x %>% str_split("_") %>% unlist
        s <- y[1]
        t <- y[2]
        pL <- parDat %>% 
            select(all_of(c(paste0("Prod_of_", t), paste0("Inh_of_", s,"To", t), 
                   paste0("Trd_of_", t, "To", s), paste0("Deg_of_", t)))) %>%
            set_names(c("p", "i", "th", "d")) %>%
            mutate(U = p, L = p*i, th = d*th) %>% 
            select(U, L, th) %>%
            mutate(k = th > L, j = th > U) %>%
            mutate(m = 2-k - j) %>% select(m) %>% unlist
        pL
    }) %>% data.frame %>% 
        set_names(str_extract(edges, "_.*") %>% str_remove("_") %>% paste0("_class")) 
    nodes <- colnames(d) %>% str_remove("_class")
    nodesOrder <- prsFile %>% filter(str_detect(Parameter, "Prod")) %>% 
        select(Parameter) %>% unlist %>% str_remove("Prod_of_")
    stateKey <- c("00", "11", "10", "01", "01-10", "10","10", "01", "01")
    names(stateKey) <- c("00","22", "20", "02", "11", "21", "10", "12", "01")
    d <- d %>% set_names(nodes) %>% select(all_of(nodesOrder)) %>% 
        set_names(paste0("N", 1:length(nodes))) %>% 
        mutate(Node = paste0(N1, N2)) %>% 
        mutate(DSGRNPhase = stateKey[Node], 
               ParamNode = Node %>% str_split("") %>% lapply(as.integer) %>% 
                   sapply(function(x) {
                        8 - 3*x[1] - x[2]
                   }),
               parIndex = parDat$parIndex) %>%
        select(parIndex, ParamNode, DSGRNPhase)
    d
}

solutionFileProcess <- function(net, tol = 1e-3, paramIDs = T, basinThresh = NULL)
{
    ### Adding parameter indices ----
    Net <- net %>% str_remove("\\d") %>% str_remove("_.*")
    nodes <- read.delim(paste0(net, ".prs")) %>% filter(str_detect(Parameter, "Prod")) %>% 
            select(Parameter) %>% unlist %>% str_remove("Prod_of_")
    solnNames <- c("parIndex", "nStatesRAC", "basin", nodes)
    solnFile <- read_delim(paste0(net, "_solution.dat"), col_names = solnNames, 
                           col_types = "d")
    bMax <- solnFile %>% filter(nStatesRAC == 1) %>% select(basin) %>% unlist %>% unique
    solnFile <- solnFile%>% mutate(basin = basin/bMax)
    pars <- solnFile$parIndex %>% unique
    
    parIds <- data.frame(parIndex = pars,
        DSGRNPhase = readLines(paste0(net, "_parameters_FPs.dat")) %>% 
            str_split("\t") %>% 
            sapply(function(x) {
                x %>% unique %>% sort %>% paste0(collapse = "-")
            }))
    


    solnFile <- merge(parIds, solnFile, by = "parIndex", all = T)
    
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
        d <- merge(sl, par, by = "parIndex", all = T) %>% 
            set_names(c("parIndex", "Node", paste0("Link", 1:i)))
        node <- 2^d$Node
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
    write_csv(states, paste0(net, "_solutionFormatter.csv"), quote = "none")
    phase <- states %>% select(parIndex, Phase) %>% distinct %>% 
        group_by(Phase) %>% summarise(Count = n(),.groups = "drop") %>%
        mutate(Frequency = Count/sum(Count))
    write_csv(phase, paste0(net, "_phases.csv"), quote = "none")
}
# sapply(c("TT", "TT1", "TT2"), solutionFileProcess)
