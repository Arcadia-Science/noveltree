# Makefile for NovelTree Docker image management
#
# This Makefile provides targets for building Docker images used in the NovelTree pipeline.
# All images are built from the repository root to ensure proper build context.

# Docker image configuration
DOCKER_PLATFORM := linux/amd64
DOCKER_ORG := arcadiascience

# Image names and versions
PHYSICOCHEMICAL_PROPS_IMAGE := $(DOCKER_ORG)/physicochemical_props
PHYSICOCHEMICAL_PROPS_TAG := 1.0.0

PHYLO_DIST_IMAGE := $(DOCKER_ORG)/phylo_dist
PHYLO_DIST_TAG := 1.0.0

# Phony targets (targets that don't represent files)
.PHONY: help docker-all docker-physicochemical-props docker-phylo-dist clean

# Default target: show help
help:
	@echo "NovelTree Docker Image Build Targets"
	@echo "====================================="
	@echo ""
	@echo "Available targets:"
	@echo "  make docker-all                  - Build all Docker images"
	@echo "  make docker-physicochemical-props - Build physicochemical properties image"
	@echo "  make docker-phylo-dist           - Build phylogenetic distance image"
	@echo "  make clean                       - Remove dangling Docker images"
	@echo ""
	@echo "Images are built for platform: $(DOCKER_PLATFORM)"
	@echo ""
	@echo "Note: Building may take 15-20 minutes depending on your system."

# Build all Docker images
docker-all: docker-physicochemical-props docker-phylo-dist
	@echo ""
	@echo "✓ All Docker images built successfully!"
	@echo ""
	@echo "Built images:"
	@echo "  - $(PHYSICOCHEMICAL_PROPS_IMAGE):$(PHYSICOCHEMICAL_PROPS_TAG)"
	@echo "  - $(PHYLO_DIST_IMAGE):$(PHYLO_DIST_TAG)"

# Build physicochemical properties Docker image
docker-physicochemical-props:
	@echo "Building physicochemical properties Docker image..."
	@echo "Image: $(PHYSICOCHEMICAL_PROPS_IMAGE):$(PHYSICOCHEMICAL_PROPS_TAG)"
	@echo "Platform: $(DOCKER_PLATFORM)"
	@echo ""
	docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(PHYSICOCHEMICAL_PROPS_IMAGE):$(PHYSICOCHEMICAL_PROPS_TAG) \
		-f docker/physicochemical_props/Dockerfile \
		.
	@echo ""
	@echo "✓ Physicochemical properties image built successfully!"

# Build phylogenetic distance Docker image
docker-phylo-dist:
	@echo "Building phylogenetic distance Docker image..."
	@echo "Image: $(PHYLO_DIST_IMAGE):$(PHYLO_DIST_TAG)"
	@echo "Platform: $(DOCKER_PLATFORM)"
	@echo ""
	docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(PHYLO_DIST_IMAGE):$(PHYLO_DIST_TAG) \
		-f docker/phylo_dist/Dockerfile \
		.
	@echo ""
	@echo "✓ Phylogenetic distance image built successfully!"

# Clean up dangling Docker images
clean:
	@echo "Removing dangling Docker images..."
	docker image prune -f
	@echo "✓ Cleanup complete!"

# Advanced: Push images to Docker Hub (requires authentication)
.PHONY: docker-push-all docker-push-physicochemical-props docker-push-phylo-dist

docker-push-physicochemical-props: docker-physicochemical-props
	@echo "Pushing $(PHYSICOCHEMICAL_PROPS_IMAGE):$(PHYSICOCHEMICAL_PROPS_TAG) to Docker Hub..."
	docker push $(PHYSICOCHEMICAL_PROPS_IMAGE):$(PHYSICOCHEMICAL_PROPS_TAG)
	@echo "✓ Image pushed successfully!"

docker-push-phylo-dist: docker-phylo-dist
	@echo "Pushing $(PHYLO_DIST_IMAGE):$(PHYLO_DIST_TAG) to Docker Hub..."
	docker push $(PHYLO_DIST_IMAGE):$(PHYLO_DIST_TAG)
	@echo "✓ Image pushed successfully!"

docker-push-all: docker-push-physicochemical-props docker-push-phylo-dist
	@echo ""
	@echo "✓ All images pushed to Docker Hub successfully!"
