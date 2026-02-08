.PHONY: all check check-fast test build document coverage install lint format clean help

all: check

help:
	@echo "prophosqua development targets:"
	@echo "  make all       - full pipeline: document -> build -> check (default)"
	@echo "  make check     - R CMD check (runs document, build first)"
	@echo "  make check-fast - R CMD check without vignettes"
	@echo "  make test      - run testthat tests (runs document first)"
	@echo "  make build     - build tarball (runs document first)"
	@echo "  make document  - generate roxygen2 docs"
	@echo "  make coverage  - code coverage report"
	@echo "  make install   - install package locally"
	@echo "  make lint      - run lintr"
	@echo "  make format    - format with air"
	@echo "  make clean     - remove build artifacts"

document:
	Rscript -e "devtools::document()"

build: document
	Rscript -e "devtools::build()"

check: build
	Rscript -e "devtools::check()"

check-fast: document
	Rscript -e "devtools::check(build_args = '--no-build-vignettes', args = '--no-vignettes')"

test: document
	Rscript -e "devtools::test()"

coverage: document
	Rscript -e "covr::package_coverage() |> print()"

install: document
	Rscript -e "devtools::install()"

lint:
	Rscript -e "lintr::lint_package()"

format:
	air format .

clean:
	rm -rf *.Rcheck
	rm -f Rplots.pdf
