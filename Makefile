.DEFAULT_GOAL := check

.PHONY: check check-fast test build build-vignettes document coverage install lint format clean help site deploy install-pre-commit-hook new_version

DOCUMENT_CMD = Rscript -e "devtools::document()"
PKG_NAME = $(shell grep '^Package:' DESCRIPTION | sed 's/Package: //')
PKG_VERSION = $(shell grep '^Version:' DESCRIPTION | sed 's/Version: //')
TARBALL = ../$(PKG_NAME)_$(PKG_VERSION).tar.gz
BUILD_CMD = Rscript -e "devtools::build(vignettes = TRUE)"
CHECK_CMD = Rscript -e "devtools::check()"
CHECK_FAST_CMD = Rscript -e "devtools::check(build_args = '--no-build-vignettes', args = '--no-vignettes')"
BUILD_VIGNETTES_CMD = Rscript -e "devtools::build_vignettes()"
TEST_CMD = Rscript -e "devtools::test()"
COVERAGE_CMD = Rscript -e "covr::package_coverage() |> print()"
INSTALL_CMD = R CMD INSTALL $(TARBALL)
LINT_CMD = Rscript -e "lintr::lint_package()"
SITE_CMD = Rscript -e "pkgdown::build_site()"
DEPLOY_CMD = Rscript -e "pkgdown::deploy_to_branch()"

help:
	@echo "prolfquapp development targets:"
	@echo "  make check     - R CMD check (runs document, build first)"
	@echo "  make check-fast - R CMD check without rebuilding vignettes during check"
	@echo "  make install   - build tarball (with vignettes) and install it"
	@echo "  make lint      - run lintr"
	@echo "  make format    - format package with air"
	@echo "  make install-pre-commit-hook - install local pre-commit hook"
	@echo "  make clean     - remove build artifacts"
	@echo "  make new_version - bump patch version, tag, and push"
	@echo ""
	@echo "Advanced:"
	@echo "  make document  - generate roxygen2 docs"
	@echo "  make build     - build tarball (with vignettes)"
	@echo "  make build-vignettes - build vignettes into inst/doc"
	@echo "  make test      - run testthat tests"
	@echo "  make coverage  - code coverage report"
	@echo "  make site      - build pkgdown site locally"
	@echo "  make deploy    - build pkgdown site and push to gh-pages"

document:
	$(DOCUMENT_CMD)

build: document
	$(BUILD_CMD)

check: build
	$(CHECK_CMD)

build-vignettes: document
	$(BUILD_VIGNETTES_CMD)
	mkdir -p inst/doc
	cp doc/*.html doc/*.Rmd doc/*.R inst/doc/ 2>/dev/null || true

check-fast: document
	$(CHECK_FAST_CMD)

test: document
	$(TEST_CMD)

coverage: document
	$(COVERAGE_CMD)

install: build
	$(INSTALL_CMD)

lint:
	$(LINT_CMD)

format:
	air format .

install-pre-commit-hook:
	cp ".githooks/pre-commit" ".git/hooks/pre-commit"
	chmod +x ".git/hooks/pre-commit"

site: document
	$(SITE_CMD)

deploy: document
	$(DEPLOY_CMD)

clean:
	rm -rf *.Rcheck
	rm -f Rplots.pdf
	rm -rf inst/doc doc Meta

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
