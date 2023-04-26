source('solutionFormat.R')
source("compareStability.R")

nList1 <- list(c(1), 1:3, 2:3, 2:4, 3:5, 4:6, 5:7, 6:8, 7:9, 8:10, 9:11, 10:12)
nList2 <- lapply(seq(2, 55, 5), function(i) {
    i:(i+9)
})
nList3 <- lapply(seq(1, 91, 10), function(i) {
    i:(i+9)
})

nList <- c(nList1, nList2, nList3)
cmd <- function(topoFile ,nL) {
    paste0("racipemt ", topoFile, " -threads 48 -num_paras 10000 -minN ", nL[1], " -maxN ", nL[length(nL)])
}

simulate <- sapply(nList, function(n) {
    net <- paste0("TT_", n[1], "To", n[length(n)])
    setwd(net)
    topoFiles <- list.files(".", ".topo")
    prs <- list.files(".", ".prs")
    param <- list.files(".", "_parameters.dat")
    soln <- list.files(".", "_solution.dat")
    nets <- topoFiles %>% str_remove(".topo")
    sapply(nets, function(Net) {
        par <- paste0(Net, "_parameters.dat")
        pr <- paste0(Net, ".prs")
        # cm <- paste0("python3 convert_racipe_parameters.py ", par, " ", pr)
        # solutionFileProcess(Net)
        compare(Net)
    })
    setwd("..")
})
