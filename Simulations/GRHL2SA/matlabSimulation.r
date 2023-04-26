source("../topo_to_matlab_DSGRNFP.R")
folders <- list.dirs(".",  full.names = F, recursive = F)
folders <- folders[str_detect(folders, "GRHL")]
scriptCommands <- c()
s <- lapply(folders, function(f) {
    setwd(f)
    nums <- f %>% str_remove(".*_") %>% str_split("To") %>% unlist %>% as.numeric
    cm <- NULL
    if (all(nums < 60)) {
        topoFiles <- list.files(".", ".topo")
        sapply(topoFiles, converter)
        cm <- c(paste0("cd ", f),
                list.files(".", "_solver.m") %>% str_remove(".m"),
                "cd ..")
    }
    setwd("..")
    return(cm)
}) %>% unlist

writeLines(s, "script.m")
