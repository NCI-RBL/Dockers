if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org")

BiocManager::install("DESeq2")

install.packages('reshape2',repos = "http://cran.us.r-project.org")
install.packages("plyr",repos = "http://cran.us.r-project.org")
