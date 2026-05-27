cores <- parallel::detectCores()
if (is.na(cores)) {
    cores <- 2L
}

cran_binary_repo <- "https://packagemanager.posit.co/cran/__linux__/focal/2021-02-01"
cran_source_repo <- "https://packagemanager.posit.co/cran/2021-02-01"
strong_deps <- c("Depends", "Imports", "LinkingTo")

options(
    repos = c(CRAN = cran_binary_repo),
    Ncpus = max(1L, cores - 1L)
)

can_load <- function(package) {
    expr <- paste0(
        "if (!requireNamespace(",
        encodeString(package, quote = "\""),
        ", quietly = TRUE)) quit(status = 1)"
    )
    status <- system2(
        file.path(R.home("bin"), "Rscript"),
        c("--vanilla", "-e", shQuote(expr)),
        stdout = FALSE,
        stderr = FALSE
    )
    identical(status, 0L)
}

install_cran <- function(packages) {
    packages <- unique(packages)
    install.packages(packages, dependencies = strong_deps)

    missing_packages <- packages[
        !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_packages) == 0) {
        return(invisible(TRUE))
    }

    options(repos = c(CRAN = cran_source_repo))
    install.packages(
        missing_packages,
        dependencies = strong_deps,
        type = "source"
    )

    still_missing <- missing_packages[!vapply(missing_packages, can_load, logical(1))]
    if (length(still_missing) > 0) {
        stop("Missing required CRAN packages: ", paste(still_missing, collapse = ", "))
    }

    options(repos = c(CRAN = cran_binary_repo))
    invisible(TRUE)
}

install_cran_source <- function(packages) {
    old_repos <- getOption("repos")
    old_ncpus <- getOption("Ncpus")
    on.exit({
        options(repos = old_repos, Ncpus = old_ncpus)
    }, add = TRUE)

    options(repos = c(CRAN = cran_source_repo), Ncpus = 1L)
    for (package in packages) {
        install.packages(
            package,
            dependencies = strong_deps,
            type = "source"
        )
        if (!can_load(package)) {
            stop("Missing required CRAN package after source install: ", package)
        }
    }

    invisible(TRUE)
}

base_cran_packages <- c(
    "BiocManager",
    "curl",
    "ggplot2",
    "matrixStats",
    "remotes",
    "stringr",
    "tidyr",
    "vroom"
)

mobster_cran_packages <- c(
    "abind",
    "ade4",
    "backports",
    "bbmle",
    "bdsmatrix",
    "boot",
    "broom",
    "car",
    "carData",
    "cellranger",
    "clipr",
    "clisymbols",
    "codetools",
    "conquer",
    "corrplot",
    "cowplot",
    "data.table",
    "doParallel",
    "entropy",
    "forcats",
    "foreign",
    "foreach",
    "fs",
    "ggforce",
    "ggplot2",
    "ggpubr",
    "ggrepel",
    "ggsci",
    "ggsignif",
    "ggthemes",
    "gridExtra",
    "GUILDS",
    "haven",
    "highr",
    "htmltools",
    "iterators",
    "knitr",
    "lattice",
    "lme4",
    "lubridate",
    "maptools",
    "markdown",
    "MASS",
    "Matrix",
    "matrixcalc",
    "MatrixModels",
    "mgcv",
    "minqa",
    "modelr",
    "mvtnorm",
    "nlme",
    "nloptr",
    "nnet",
    "numDeriv",
    "openxlsx",
    "pbkrtest",
    "pixmap",
    "plyr",
    "poilog",
    "polyclip",
    "polynom",
    "pracma",
    "prettydoc",
    "quantreg",
    "RcppArmadillo",
    "RcppEigen",
    "readr",
    "readxl",
    "rematch",
    "reprex",
    "reshape2",
    "rio",
    "rmarkdown",
    "rstatix",
    "rvest",
    "sads",
    "segmented",
    "seqinr",
    "sp",
    "SparseM",
    "statmod",
    "tidyverse",
    "tinytex",
    "tweenr",
    "VGAM",
    "viridis",
    "xfun",
    "yaml",
    "zip"
)

install_cran(c(base_cran_packages, mobster_cran_packages))
install_cran_source(c("igraph", "graphlayouts", "tidygraph", "ggraph"))

if (utils::packageVersion("matrixStats") <= "0.57.0") {
    options(repos = c(CRAN = cran_source_repo))
    install.packages("matrixStats", dependencies = strong_deps, type = "source")
    options(repos = c(CRAN = cran_binary_repo))
}
if (utils::packageVersion("matrixStats") <= "0.57.0") {
    stop("matrixStats must be newer than 0.57.0 for Bioconductor 3.12")
}

BiocManager::install(version = "3.12", ask = FALSE, update = FALSE)
BiocManager::install(
    c("XVector", "GenomicRanges", "VariantAnnotation"),
    ask = FALSE,
    update = FALSE
)

github_packages <- c(
    "caravagn/pio@2cfe575",
    "caravagn/easypar@226e216",
    "caravagn/ctree@df8dc83",
    "im3sanger/dndscv@69007c2",
    "caravagnalab/mobster@85c898f"
)
for (github_package in github_packages) {
    remotes::install_github(
        github_package,
        dependencies = FALSE,
        upgrade = "never"
    )
}

required_packages <- c("XVector", "GenomicRanges", "VariantAnnotation", "mobster")
missing_packages <- required_packages[
    !vapply(required_packages, can_load, logical(1))
]
if (length(missing_packages) > 0) {
    stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}
