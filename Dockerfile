FROM r-base:4.4.0

RUN apt-get update \
  && apt-get install -y pkg-config bzip2 curl openssl \
  libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
  libfontconfig1-dev libfreetype-dev libexpat1 libexpat1-dev \
  cmake libharfbuzz-dev libfribidi-dev libxml2-dev libtiff-dev \
  && rm -rf /var/lib/apt/lists/*

# ENV PKG_CONFIG_PATH=/usr/lib/x86_64-linux-gnu/pkgconfig
RUN R -e "install.packages('devtools', Ncpus=4)"
RUN R -e "library(devtools);install_github('PheWAS/PheWAS')"

