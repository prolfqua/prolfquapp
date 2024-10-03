ARG R_VERSION=4.4.1

FROM r-base:${R_VERSION} AS base
ARG TARGETPLATFORM
ARG QUARTO_VERSION=1.5.57

RUN apt-get update \
  && apt-get install -y pandoc gdebi \
  && rm -rf /var/lib/apt/lists/*

RUN apt-get update

RUN if [ "$TARGETPLATFORM" = "linux/arm64" ]; then ARCHITECTURE=arm64; else ARCHITECTURE=amd64; fi \
  && wget "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${ARCHITECTURE}.deb" -O /tmp/quarto.deb \
  && gdebi --non-interactive /tmp/quarto.deb \
  && rm /tmp/quarto.deb



FROM base AS build
SHELL ["/bin/bash", "-c"]

RUN apt-get update \
  && apt-get upgrade -y \
  && apt-get install -y libcurl4-openssl-dev cmake libglpk-dev libxml2-dev libfontconfig1-dev libfreetype6-dev \
  && rm -rf /var/lib/apt/lists/*
ENV R_LIBS_USER=/opt/r-libs-site
RUN mkdir -p /opt/r-libs-site

RUN R -e 'options(warn=2); install.packages("pak", repos = "https://stat.ethz.ch/CRAN/")'
RUN R -e 'options(warn=2); pak::pkg_install(c("any::seqinr", "any::prozor", "any::logger", "any::lubridate", "git::https://gitlab.bfabric.org/wolski/prolfquadata.git", "github::fgcz/prolfqua"))'
COPY ./DESCRIPTION /opt/prolfqua/DESCRIPTION
RUN R -e 'options(warn=2); pak::local_install_deps("/opt/prolfqua", upgrade = FALSE)'
COPY . /opt/prolfqua
RUN R -e 'options(warn=2); pak::pkg_install("/opt/prolfqua", upgrade = FALSE)'



FROM base
COPY --from=build /opt/r-libs-site /opt/r-libs-site
RUN mkdir -p /tmp/quarto-cache && chmod 0777 /tmp/quarto-cache
ENV XDG_CACHE_HOME=/tmp/quarto-cache
ENV R_LIBS_USER=/opt/r-libs-site
ENV PATH="/opt/r-libs-site/prolfquapp/application/bin:/root/.local/bin:${PATH}"
ENTRYPOINT ["/bin/bash"]
