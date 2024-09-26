ARG R_VERSION=4.4.1

FROM r-base:${R_VERSION} AS base
RUN apt-get update \
  && apt-get install -y pandoc \
  && rm -rf /var/lib/apt/lists/*



FROM base AS build
SHELL ["/bin/bash", "-c"]

RUN apt-get update \
  && apt-get install -y libcurl4-openssl-dev cmake libglpk-dev libxml2-dev libfontconfig1-dev libfreetype6-dev \
  && rm -rf /var/lib/apt/lists/*
ENV R_LIBS_USER=/opt/r-libs-site
RUN mkdir -p /opt/r-libs-site

RUN R -e 'options(warn=2); install.packages("pak", repos = "https://stat.ethz.ch/CRAN/")'
RUN R -e 'options(warn=2); pak::pkg_install(c("any::seqinr", "any::prozor", "any::logger", "git::https://gitlab.bfabric.org/wolski/prolfquadata.git", "github::fgcz/prolfqua"))'
COPY ./DESCRIPTION /opt/prolfqua/DESCRIPTION
RUN R -e 'options(warn=2); pak::local_install_deps("/opt/prolfqua", upgrade = FALSE)'
COPY . /opt/prolfqua
RUN R -e 'options(warn=2); pak::pkg_install("/opt/prolfqua", upgrade = FALSE)'



FROM base
COPY --from=build /opt/r-libs-site /opt/r-libs-site
ENV R_LIBS_USER=/opt/r-libs-site
ENV PATH="/opt/r-libs-site/prolfquapp/application/bin:${PATH}"
ENTRYPOINT ["/bin/bash"]
