# syntax=docker/dockerfile:1
ARG R_VERSION=4.6.1
ARG URGENT_BASE_IMAGE=ghcr.io/prolfqua/prolfquapp:latest

# Pin the multi-architecture Posit Ubuntu Noble image by digest so the build and
# runtime stages always share the same R version and operating-system ABI.
FROM posit/r-base:${R_VERSION}-noble@sha256:632f7e3e493c10817a17b52065efb9f48a977ecc1f470e236323f36f4ed1db19 AS base
ARG TARGETPLATFORM
ARG TARGETARCH
ARG QUARTO_VERSION=1.5.57
ARG RSPM_R_VERSION=4.6

# Use Posit's native Ubuntu Noble binaries for both target architectures.
RUN case "$TARGETARCH" in \
    amd64) RSPM_ARCH=x86_64 ;; \
    arm64) RSPM_ARCH=aarch64 ;; \
    *) echo "Unsupported target architecture for RSPM: $TARGETARCH" >&2; exit 1 ;; \
  esac \
  && RSPM_REPO="https://packagemanager.posit.co/cran/latest/bin/linux/noble-${RSPM_ARCH}/${RSPM_R_VERSION}" \
  && printf 'options(repos = c(CRAN = "%s"))\n' "$RSPM_REPO" >> "$(R RHOME)/etc/Rprofile.site"

RUN if [ "$TARGETPLATFORM" = "linux/arm64" ]; then ARCHITECTURE=arm64; else ARCHITECTURE=amd64; fi \
  && apt-get update \
  && apt-get install -y --no-install-recommends ca-certificates pandoc wget \
  && wget "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${ARCHITECTURE}.deb" -O /tmp/quarto.deb \
  && apt-get install -y --no-install-recommends /tmp/quarto.deb \
  && rm /tmp/quarto.deb \
  && rm -rf /var/lib/apt/lists/*



FROM base AS build
SHELL ["/bin/bash", "-c"]

RUN apt-get update \
  && apt-get install -y libcurl4-openssl-dev cmake libglpk-dev libxml2-dev libfontconfig1-dev libfreetype6-dev \
  && rm -rf /var/lib/apt/lists/*
ENV R_LIBS_USER=/opt/r-libs-site
RUN mkdir -p /opt/r-libs-site

# Install arrow with its full feature set (compression codecs incl. zstd) before
# dependency resolution. RSPM normally supplies a native Noble binary; the build
# flags preserve full codec support if it falls back to a source installation.
ENV LIBARROW_MINIMAL=false
ENV ARROW_WITH_ZSTD=ON
RUN R -e 'options(warn=2); install.packages("arrow"); stopifnot(arrow::arrow_info()$capabilities[["zstd"]])'

RUN R -e 'options(warn=2); install.packages("pak")'
COPY <<'EOF' /tmp/install_cran_binary_deps.R
args <- commandArgs(trailingOnly = TRUE)
mode <- args[[1]]
targets <- args[-1]
deps <- switch(
  mode,
  refs = pak::pkg_deps(targets, upgrade = FALSE),
  local = pak::local_deps(targets[[1]], upgrade = FALSE),
  description = {
    description <- read.dcf(targets[[1]], fields = "Remotes")
    remotes <- trimws(strsplit(description[1, "Remotes"], ",", fixed = TRUE)[[1]])
    stopifnot(length(remotes) > 0L, all(nzchar(remotes)))
    pak::pkg_deps(remotes, upgrade = FALSE)
  },
  stop("Unknown dependency mode: ", mode)
)
available <- rownames(available.packages())
installed <- rownames(installed.packages())
packages <- sort(setdiff(intersect(unique(deps$package), available), installed))
if (length(packages) > 0L) {
  message("Installing ", length(packages), " resolved CRAN dependencies as native binaries")
  install.packages(packages)
}
EOF
RUN Rscript /tmp/install_cran_binary_deps.R refs \
  any::seqinr any::prozor any::logger any::lubridate \
  github::fgcz/prolfqua github::prolfqua/prolfquasaint
RUN R -e 'options(warn=2); pak::pkg_install(c("any::seqinr", "any::prozor", "any::logger", "any::lubridate", "github::fgcz/prolfqua", "github::prolfqua/prolfquasaint"))'
COPY ./DESCRIPTION /opt/prolfqua/DESCRIPTION
RUN Rscript /tmp/install_cran_binary_deps.R local /opt/prolfqua
RUN R -e 'options(warn=2); pak::local_install_deps("/opt/prolfqua", upgrade = FALSE)'
COPY . /opt/prolfqua
RUN R -e 'options(warn=2); install.packages(c("knitr", "rmarkdown", "DT", "gridExtra", "KernSmooth", "plotly", "Rfit"))'
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
stopifnot(file.exists(file.path(out, "Grp2Analysis_V2_SE_tabset.html")))
cat("Quarto SE report rendered OK ->", file.path(out, "Grp2Analysis_V2_SE_tabset.html"), "\n")
EOF

RUN Rscript /tmp/check_vignettes.R \
  && Rscript /tmp/check_quarto.R

# Validate data.table loads correctly at build time
RUN Rscript -e "cat('Testing data.table load...\\n'); library(data.table); cat('data.table loaded successfully.\\n')"


FROM base AS full
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
COPY --from=build /tmp/check_vignettes.R /tmp/check_quarto.R /tmp/install_cran_binary_deps.R /opt/checks/
RUN mkdir -p /tmp/quarto-cache && chmod 0777 /tmp/quarto-cache
ENV XDG_CACHE_HOME=/tmp/quarto-cache
ENV R_LIBS_USER=/opt/r-libs-site
RUN for dir in /opt/.TinyTeX/bin/*/; do ln -sf $dir* /usr/local/bin/; done
ENV PATH="/opt/r-libs-site/prolfquapp/application/bin:/root/.local/bin:${PATH}"
ENTRYPOINT ["/bin/bash"]


# Urgent patch releases reuse the environment from their matching X.Y.0 full
# image, refresh every package declared in Remotes, and rebuild prolfquapp. The
# workflow verifies that Dockerfile and package dependencies have not changed
# since that full-release tag.
FROM ${URGENT_BASE_IMAGE} AS urgent-dependencies
COPY DESCRIPTION /tmp/prolfquapp-DESCRIPTION
RUN Rscript /opt/checks/install_cran_binary_deps.R description /tmp/prolfquapp-DESCRIPTION \
  && R -e 'description <- read.dcf("/tmp/prolfquapp-DESCRIPTION", fields = "Remotes"); remotes <- trimws(strsplit(description[1, "Remotes"], ",", fixed = TRUE)[[1]]); stopifnot(length(remotes) > 0L, all(nzchar(remotes))); message("Updating urgent release remotes: ", paste(remotes, collapse = ", ")); pak::pkg_install(remotes, upgrade = FALSE)'


FROM urgent-dependencies AS urgent-build
ARG RELEASE_VERSION
COPY . /opt/prolfqua
RUN set -eu; \
  actual_version="$(sed -n 's/^Version:[[:space:]]*//p' /opt/prolfqua/DESCRIPTION)"; \
  if [ "$actual_version" != "$RELEASE_VERSION" ]; then \
    echo "DESCRIPTION version $actual_version does not match release tag $RELEASE_VERSION" >&2; \
    exit 1; \
  fi; \
  cd /tmp; \
  R CMD build /opt/prolfqua --no-manual; \
  set -- /tmp/prolfquapp_*.tar.gz; \
  if [ "$#" -ne 1 ]; then \
    echo "Expected exactly one prolfquapp source tarball, found $#" >&2; \
    exit 1; \
  fi; \
  mv "$1" /tmp/prolfquapp-release.tar.gz


FROM urgent-dependencies AS urgent
ARG RELEASE_VERSION
COPY --from=urgent-build /tmp/prolfquapp-release.tar.gz /tmp/prolfquapp-release.tar.gz
RUN set -eu; \
  R CMD INSTALL /tmp/prolfquapp-release.tar.gz; \
  rm /tmp/prolfquapp-release.tar.gz; \
  installed_version="$(Rscript -e 'library(prolfquapp); cat(as.character(packageVersion("prolfquapp")))')"; \
  if [ "$installed_version" != "$RELEASE_VERSION" ]; then \
    echo "Installed prolfquapp version $installed_version does not match release tag $RELEASE_VERSION" >&2; \
    exit 1; \
  fi
