
FROM rocker/r-ver:4.3.0

# System dependencies for R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('optparse', repos='https://cloud.r-project.org/')"

RUN mkdir -p /opt/scripts

COPY scripts/ /opt/scripts/

RUN chmod +x /opt/scripts/*.r

ENV PATH="/opt/scripts:${PATH}"