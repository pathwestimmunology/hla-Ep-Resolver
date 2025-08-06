# Get base image
FROM rocker/r-base:4.5.0

# Set environment variables
ENV TZ="Australia/Perth"
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# For support, contact:
LABEL maintainer="Fredrick Mobegi <fredrick.mobegi@health.wa.gov.au>"

# system libraries of general use
## install debian packages
# Install libraries that are required by Shiny and related R packages.
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libxml2-dev \
    libcurl4-openssl-dev && \
    apt-get clean

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    apt-utils \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    sudo \
    libxml2 \
    aptitude \
    libcurl4-gnutls-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libxml2-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev &&  \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean

## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean

# Set the CRAN mirror
ARG CRAN_MIRROR=https://cran.rstudio.com/

# Set the list of R packages to install
ARG R_PACKAGES="\
    pacman \
    pak \
    vroom \
    plyr \
    RCurl \
    remotes \
    BiocManager \
    shiny \
    shinyjs \
    shinyBS \
    bslib \
    shinyalert \
    shinythemes \
    shinydashboard \
    shinyalert \
    shinycssloaders \
    DataExplorer \
    modeldata \
    DT \
    dplyr \
    plyr \
    MASS \
    stringr \
    data.table \
    foreach \
    openxlsx \
    zoo \
    doParallel"

# Install R packages
RUN R -e "chooseCRANmirror(graphics=FALSE, ind=6)"

RUN for pkg in $R_PACKAGES; do \
    R -e "install.packages('$pkg', repos='$CRAN_MIRROR')"; \
    if [ $? -ne 0 ]; then \
    echo "Error: Failed to install package '$pkg'."; \
    exit 1; \
    fi; \
    done

RUN R -e "pak::pak('tidyverse/tidyverse')"

RUN export R_LIBS="/usr/local/lib/R/site-library"
ENV R_LIBS="/usr/local/lib/R/site-library"
ENV R_LIBS_SITE="/usr/local/lib/R/site-library"

# Install XML
# RUN R -e 'remotes::install_version("XML", version = "3.99-0.3")'

# Install 'InterMineR' package || requires RCurl package to work properly
RUN wget https://cran.r-project.org/src/contrib/RCurl_1.98-1.16.tar.gz && \
    tar xvzf RCurl_1.98-1.16.tar.gz && \
    R CMD INSTALL RCurl

RUN R -e "remotes::install_github('https://github.com/intermine/InterMineR')"

# Install 'immunotation' from bioconductor
RUN R -e "BiocManager::install(ask = F)" && R -e "BiocManager::install('immunotation', ask = F)"

# Set the options in the Rprofile.site file which are going to be loaded by the R session.
# These options specify the Shiny host and port that runApp will use.
RUN echo "options(shiny.port = 3839, shiny.host = '0.0.0.0', shiny.launch.browser = FALSE)" > /usr/lib/R/etc/Rprofile.site

RUN apt-get clean && \
    apt-get autoremove -y

# Create a Linux group and user, both called epletmatch.
# This user will have access to the epletmatch instead of the default root user:
# Very helpful in reducing vulnerabilities of the containerized applications.

## Doesn't work on new Ubuntu version.. updated below
# RUN addgroup --system epletmatch \
#     && adduser --system --ingroup epletmatch epletmatch
RUN groupadd --system epletmatch \
    && useradd --system --gid epletmatch epletmatch

# Set the working directory for subsequent instructions
WORKDIR /home/epletmatch

# Copy new files or directories from the source of our Shiny app and adds them to the file system of the container at the destination
# path (. refers to the current work directory defined at the previous step):
COPY app/ .

RUN mkdir -p /home/epletmatch/data && \
    mkdir -p /home/epletmatch/www

VOLUME /home/epletmatch/data
VOLUME /home/epletmatch/www

COPY data ./data
COPY www ./www
RUN ls -la .

# Set permissions for the epletmatch user:
RUN chown epletmatch:epletmatch -R /home/epletmatch

# Set the user name (or UID) and optionally the user group (or GID) to use when running the image:
USER epletmatch

# Tell Docker which ports the container listens on at runtime.
# Set this to the Shiny port defined in the Rprofile.site file above:
EXPOSE 3839

# Specify the executables. Only one CMD is allowed in a Dockerfile.
# If more than 1, only the last one will take effect.
CMD ["R", "-e", "shiny::runApp('/home/epletmatch')"]
