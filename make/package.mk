SITE_CMD = Rscript -e "altdoc::render_docs()"

DOCKER_IMAGE ?= prolfqua/$(PKG_NAME)
DOCKER_TAG ?= dev

.PHONY: sync-quarto-assets docker-build docker-check

build build-vignettes test check-fast vignette precommit-prepare: sync-quarto-assets

sync-quarto-assets:
	Rscript data-raw/sync_quarto_assets.R

docker-build: sync-quarto-assets
	docker build -t $(DOCKER_IMAGE):$(DOCKER_TAG) -f Dockerfile .

docker-check:
	docker run --rm $(DOCKER_IMAGE):$(DOCKER_TAG) -c \
	  'Rscript /opt/checks/check_vignettes.R && Rscript /opt/checks/check_quarto.R'

help-package:
	@echo ""
	@echo "Package-specific:"
	@echo "  make sync-quarto-assets - synchronize FGCZ Quarto assets"
	@echo "  make docker-build       - build $(DOCKER_IMAGE):$(DOCKER_TAG)"
	@echo "  make docker-check       - run vignette and Quarto checks in the image"
