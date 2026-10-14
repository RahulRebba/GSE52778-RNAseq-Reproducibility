# Prototype container definition for the GSE52778 project.
# The validated workflow for this repository is the local
# Python + R/renv + SQLite workflow documented in README.md.
# This container has not been validated as the final production runtime.

FROM bioconductor/bioconductor_docker:RELEASE_3_21

RUN R -q -e 'BiocManager::install(c("airway","DESeq2","clusterProfiler","org.Hs.eg.db","enrichplot"), ask=FALSE, update=FALSE); install.packages(c("ggplot2","pheatmap"))'

RUN apt-get update && \
    apt-get install -y python3 python3-pip sqlite3 && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install --break-system-packages pandas

WORKDIR /project
COPY . /project

CMD ["bash"]