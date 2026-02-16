.PHONY: all check check-fast test build document coverage install lint format clean help new_version pkgdown

all: check

help:
	@echo "prolfquapp development targets:"
	@echo "  make all         - full pipeline: document -> build -> check (default)"
	@echo "  make check       - R CMD check (runs document, build first)"
	@echo "  make check-fast  - R CMD check without vignettes"
	@echo "  make test        - run testthat tests (runs document first)"
	@echo "  make build       - build tarball (runs document first)"
	@echo "  make document    - generate roxygen2 docs"
	@echo "  make coverage    - code coverage report"
	@echo "  make install     - install package locally"
	@echo "  make lint        - run lintr"
	@echo "  make format      - format with air"
	@echo "  make clean       - remove build artifacts"
	@echo "  make new_version - bump patch version, tag, and push"
	@echo "  make pkgdown     - build pkgdown site locally"

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

new_version:
	@CURRENT=$$(grep '^Version:' DESCRIPTION | sed 's/Version: //'); \
	MAJOR=$$(echo $$CURRENT | cut -d. -f1); \
	MINOR=$$(echo $$CURRENT | cut -d. -f2); \
	PATCH=$$(echo $$CURRENT | cut -d. -f3); \
	NEW_PATCH=$$((PATCH + 1)); \
	NEW_VERSION="$$MAJOR.$$MINOR.$$NEW_PATCH"; \
	echo "Bumping version: $$CURRENT -> $$NEW_VERSION"; \
	sed -i '' "s/^Version: .*/Version: $$NEW_VERSION/" DESCRIPTION; \
	git add DESCRIPTION; \
	git commit -m "new version $$NEW_VERSION"; \
	git tag "$$NEW_VERSION"; \
	git push && git push --tags; \
	echo "Released $$NEW_VERSION"

pkgdown:
	Rscript -e "pkgdown::build_site()"
