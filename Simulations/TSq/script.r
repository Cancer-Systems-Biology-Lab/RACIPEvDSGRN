source('solutionFormat.R')
source("compareStability.R")

nList1 <- list(c(1), 1:3, 2:3, 2:4, 3:5, 4:6, 5:7, 6:8, 7:9, 8:10, 9:11, 10:12)
nList2 <- lapply(seq(2, 55, 5), function(i) {
    i:(i+9)
})
nList <- lapply(seq(1, 91, 10), function(i) {
    i:(i+9)
})
nList <- c(nList1, nList2, nList)

# nList <- c(nList1, nList2, nList3)
cmd <- function(topoFile ,nL) {
    paste0("racipemt ", topoFile, " -threads 48 -num_paras 10000 -minN ", nL[1], " -maxN ", nL[length(nL)])
}

simulate <- sapply(nList, function(n) {
    net <- paste0("TSq_", n[1], "To", n[length(n)])
    net1 <- paste0("TSq1_", n[1], "To", n[length(n)])
    net2 <- paste0("TSq2_", n[1], "To", n[length(n)])
    dir.create(net, showWarnings = F)
    file.copy("TSq.topo", paste0(net, "/", net, ".topo"))
    file.copy(paste0(net, "/", net, ".topo"), paste0(net, "/", net2, ".topo"))
    file.copy(paste0(net, "/", net, ".topo"), paste0(net, "/", net1, ".topo"))
    setwd(net)
    topoFiles <- list.files(".", ".topo")
    commands <- sapply(topoFiles, function(topoFile) {
        cmd(topoFile, n)
    })
    for (command in commands) {
        system(command)
    }
    prs <- list.files(".", ".prs")
    param <- list.files(".", "_parameters.dat")
    soln <- list.files(".", "_solution.dat")
    keep <- c(topoFiles, prs, param, soln)
    filsAll <- list.files(".")
    filsDel <- filsAll[!(filsAll %in% keep)]
    sapply(filsDel, file.remove)
    nets <- topoFiles %>% str_remove(".topo")
    sapply(nets, function(Net) {
        par <- paste0(Net, "_parameters.dat")
        pr <- paste0(Net, ".prs")
        cm <- paste0("python3 ../convert_racipe_parameters.py ", par, " ", pr)
        system(cm)
        solutionFileProcess(Net)
        compare(Net)
    })
    setwd("..")
})
