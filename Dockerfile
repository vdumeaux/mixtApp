# Use your updated base image
FROM kvik-r-updated

# Update package lists
RUN apt-get update

# Install necessary packages
RUN apt-get update && apt-get install -y \
    default-jre \
    default-jdk \
    libmariadb-dev \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    git \
    wget

# Install the 'pak' package
RUN R -e "install.packages('pak', repos = 'https://cloud.r-project.org')"

RUN apt-get update && DEBIAN_FRONTEND=noninteractive bash -c "apt-get install -y $$(Rscript -e \"cat(pak::pkg_system_requirements(c('systemfonts', 'textshaping', 'ragg', 'pkgdown', 'devtools'), os = 'debian', os_release = 'bullseye'))\") && rm -rf /var/lib/apt/lists/*"

# Install R packages from CRAN, including missing dependencies
RUN R -e 'install.packages("pkgdown", repos="https://cloud.r-project.org")'
# RUN R -e 'install.packages(c("ggplot2", "Hmisc", "Rcpp", "roxygen2", "jsonlite", "igraph", "dplyr", "parallel", "colorspace", "ic10", "network", "GGally", "sna", "animation", "latticeExtra", "reshape"), repos="https://cloud.r-project.org")'

# Install Bioconductor packages, including missing dependencies
# RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("genefu", "WGCNA", "impute", "preprocessCore", "GO.db", "illuminaHumanv4.db", "illuminaHumanv3.db", "hgu133a.db", "breastCancerVDX"), ask=FALSE)'

# Install the 'mixtR' package from GitHub
# RUN R -e 'devtools::install_github("vdumeaux/mixtR")'

# Copy your application code into the container
# COPY . /mixtApp
# RUN R CMD INSTALL /mixtApp

# Set the working directory (adjust this if needed)
# WORKDIR /go/src/github.com/fjukstad/kvik/r/examples

# (Optional) Expose necessary ports or set other configurations

