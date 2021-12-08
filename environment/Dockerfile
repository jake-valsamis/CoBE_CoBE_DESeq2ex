# hash:sha256:436478a3608ddc2ddf7d599b2fb4d6d7cef0ac42f8b6daf388dc8e1aee4b7fd0
FROM registry.codeocean.com/codeocean/r-studio:1.4.1106-r4.0.5-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "DESeq2", \
        "airway", \
        "apeglm", \
        "ashr", \
        "knitr", \
        "pasilla", \
        "pheatmap", \
        "readr", \
        "tximeta", \
        "tximport", \
        "tximportData", \
        "vsn" \
    ))' # Original versions: 1.34.0 1.13.0 1.16.0 2.2-47 1.36 1.21.0 1.0.12 2.0.2 1.12.0 1.22.0 1.21.0 3.62.0
