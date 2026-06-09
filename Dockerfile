# syntax=docker/dockerfile:1
ARG R_VERSION=4.5.2

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
RUN R -e 'options(warn=2); pak::pkg_install(c("any::seqinr", "any::prozor", "any::logger", "any::lubridate", "github::fgcz/prolfqua", "github::prolfqua/saintexpress", "github::prolfqua/prolfquasaint"))'
COPY ./DESCRIPTION /opt/prolfqua/DESCRIPTION
RUN R -e 'options(warn=2); pak::local_install_deps("/opt/prolfqua", upgrade = FALSE)'
COPY . /opt/prolfqua
RUN R -e 'options(warn=2); install.packages(c("knitr", "rmarkdown", "DT", "gridExtra", "KernSmooth", "plotly"), repos = "https://stat.ethz.ch/CRAN/")'
RUN cd /tmp \
  && R CMD build /opt/prolfqua --no-manual \
  && R CMD INSTALL prolfquapp_*.tar.gz
# Write the build-time check scripts to /tmp using BuildKit heredoc COPY,
# so multi-line R code is preserved verbatim (no shell quoting hazards).
COPY <<'EOF' /tmp/check_vignettes.R
# Re-render every Rmd vignette with rmarkdown::render() to prove that the
# installed image has every package + tool each vignette needs.
doc <- system.file("doc", package = "prolfquapp", mustWork = TRUE)
rmds <- list.files(doc, pattern = "\\.Rmd$", full.names = TRUE)
if (length(rmds) == 0) stop("No Rmd vignettes installed under doc/")
out_dir <- tempfile("vignette_check_")
dir.create(out_dir, recursive = TRUE)
for (rmd in rmds) {
  cat("Rendering: ", basename(rmd), "\n", sep = "")
  out <- rmarkdown::render(rmd, output_dir = out_dir, quiet = TRUE, envir = new.env())
  stopifnot(file.exists(out))
}
cat("All ", length(rmds), " Rmd vignettes rendered via rmarkdown::render().\n", sep = "")
EOF

COPY <<'EOF' /tmp/check_quarto.R
# Render the Quarto SE report template against the bundled SE fixture to
# verify the quarto CLI + every package the template's setup chunk needs.
out <- tempfile("qmd_check_")
dir.create(out, recursive = TRUE)
se <- system.file("extdata", "3106962.rds", package = "prolfquapp", mustWork = TRUE)
prolfquapp:::render_quarto_se_report(se_file = se, output_dir = out)
stopifnot(file.exists(file.path(out, "Grp2Analysis_V2_SE.html")))
cat("Quarto SE report rendered OK ->", file.path(out, "Grp2Analysis_V2_SE.html"), "\n")
EOF

RUN Rscript /tmp/check_vignettes.R \
  && Rscript /tmp/check_quarto.R

# Validate data.table loads correctly at build time
RUN Rscript -e "cat('Testing data.table load...\\n'); library(data.table); cat('data.table loaded successfully.\\n')"


FROM base
ARG TARGETPLATFORM
ENV HOME=/home/user
ARG TINYTEX_VERSION=2026.05
RUN if [ "$TARGETPLATFORM" = "linux/arm64" ]; then ARCH=linux-arm64; else ARCH=linux-x86_64; fi \
  && wget "https://github.com/rstudio/tinytex-releases/releases/download/v${TINYTEX_VERSION}/TinyTeX-1-${ARCH}-v${TINYTEX_VERSION}.tar.xz" -O /tmp/tinytex.tar.xz \
  && tar -xJf /tmp/tinytex.tar.xz -C /opt \
  && rm /tmp/tinytex.tar.xz
RUN mkdir -p /home/user && chmod -R 777 /home/user
COPY --from=build /opt/r-libs-site /opt/r-libs-site
# Bake the build-time check scripts into the deploy image so
# `make docker-check` can re-run them against a published image without
# needing a rebuild.
COPY --from=build /tmp/check_vignettes.R /tmp/check_quarto.R /opt/checks/
RUN mkdir -p /tmp/quarto-cache && chmod 0777 /tmp/quarto-cache
ENV XDG_CACHE_HOME=/tmp/quarto-cache
ENV R_LIBS_USER=/opt/r-libs-site
RUN for dir in /opt/.TinyTeX/bin/*/; do ln -sf $dir* /usr/local/bin/; done
ENV PATH="/opt/r-libs-site/prolfquapp/application/bin:/root/.local/bin:${PATH}"
ENTRYPOINT ["/bin/bash"]
