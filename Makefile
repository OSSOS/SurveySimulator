
REPO = images.canfar.net
PROJECT = uvickbos
DEVNAME = ssim
VERSION := $(shell python -m setuptools_scm --strip-dev)

# CANFAR astroml base tag the unified image is built on (override: make ASTROML_TAG=26.04 ...).
ASTROML_TAG ?= 26.06

NAME = $(REPO)/$(PROJECT)/$(DEVNAME)

# Build the single unified CANFAR/skaha + Cursor image. One image runs in all three
# skaha session types (notebook/desktop-app/headless) and is the Cursor/devcontainer base.
build: Dockerfile
	docker build --build-arg ASTROML_TAG=$(ASTROML_TAG) -t $(NAME):$(VERSION) -t $(NAME):latest -f Dockerfile .

# Backwards-compatible alias for the previous `production` target.
production: build

# Push the built image to the CANFAR Harbor registry (no rebuild).
# Run `docker login images.canfar.net` first.
push:
	docker push $(NAME):$(VERSION)
	docker push $(NAME):latest

# Build then push a release to canfar.net.
deploy: build push

# Build then print how to run the image locally as a desktop-app xterm.
dev: build
	echo "docker run --rm -it $(NAME):$(VERSION) xterm"

.PHONY: build production push deploy dev clean
clean:
	\rm -rf build
