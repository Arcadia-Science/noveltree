# Makefile for NovelTree Docker Image Management
#
# This Makefile provides targets for building Docker images used in the NovelTree pipeline.
# Images that reference files outside their docker/ subdirectory are built from the repository root.
# All other images are built from their respective docker/ subdirectories for cleaner, faster builds.

# Docker configuration
DOCKER_PLATFORM := linux/amd64
DOCKER_ORG := arcadiascience

# =============================================================================
# Image names and versions
# =============================================================================

# Images requiring root context (COPY files from outside their docker/ directory)
PROTEIN_PROPERTIES_IMAGE := $(DOCKER_ORG)/protein_properties
PROTEIN_PROPERTIES_TAG := 1.0.0
ZOOGLE_IMAGE := $(DOCKER_ORG)/zoogle
ZOOGLE_TAG := 1.0.0

# Images built from their docker/ subdirectory
ASTEROID_IMAGE := $(DOCKER_ORG)/asteroid_3aae117d-disco_20e10c33
ASTEROID_TAG := 1.0.0
BIOSERVICES_IMAGE := $(DOCKER_ORG)/bioservices_1.10.0
BIOSERVICES_TAG := 1.0.0
CIALIGN_IMAGE := $(DOCKER_ORG)/cialign_1.1.0
CIALIGN_TAG := 1.0.0
CLIPKIT_IMAGE := $(DOCKER_ORG)/clipkit_2.1.1-seqmagick_0.8.4
CLIPKIT_TAG := 1.0.0
COGEQC_IMAGE := $(DOCKER_ORG)/cogeqc_1.2.1
COGEQC_TAG := 1.0.0
FAMSA_IMAGE := $(DOCKER_ORG)/famsa_2.0.0
FAMSA_TAG := 1.0.0
FASTTREE_IMAGE := $(DOCKER_ORG)/fasttree_2.1.11
FASTTREE_TAG := 1.0.0
GENERAX_IMAGE := $(DOCKER_ORG)/generax_56f3ed0
GENERAX_TAG := 1.1.3
IQTREE_IMAGE := $(DOCKER_ORG)/iqtree_2.2.0.5
IQTREE_TAG := 1.0.0
ORTHOFINDER_IMAGE := $(DOCKER_ORG)/orthofinder_2.5.4
ORTHOFINDER_TAG := 1.0.0
PHYLO_PROFILES_IMAGE := $(DOCKER_ORG)/phylo_profiles
PHYLO_PROFILES_TAG := 1.0.0
RBASE_IMAGE := $(DOCKER_ORG)/rbase_4.2.2
RBASE_TAG := 1.0.0
SELECT_INFLATION_IMAGE := $(DOCKER_ORG)/select_mcl_inflation_params_08302023
SELECT_INFLATION_TAG := 1.0.0
WITCH_IMAGE := $(DOCKER_ORG)/witch_0.3.0
WITCH_TAG := 1.0.0
PREPROCESS_PROTEOMES_IMAGE := $(DOCKER_ORG)/preprocess_proteomes
PREPROCESS_PROTEOMES_TAG := 1.1.0
PREQUAL_IMAGE := $(DOCKER_ORG)/prequal
PREQUAL_TAG := 1.0.0

# =============================================================================
# Phony targets
# =============================================================================

.PHONY: help docker-all \
	docker-protein-properties docker-zoogle \
	docker-asteroid docker-bioservices docker-cialign docker-clipkit \
	docker-cogeqc docker-famsa docker-fasttree docker-generax \
	docker-iqtree docker-orthofinder \
	docker-phylo-profiles docker-rbase docker-select-inflation docker-witch \
	docker-preprocess-proteomes docker-prequal \
	push-all push-protein-properties push-zoogle \
	push-asteroid push-bioservices push-cialign push-clipkit \
	push-cogeqc push-famsa push-fasttree push-generax \
	push-iqtree push-orthofinder \
	push-phylo-profiles push-rbase push-select-inflation push-witch \
	push-preprocess-proteomes push-prequal \
	clean

# =============================================================================
# Help target
# =============================================================================

help:
	@echo "NovelTree Docker Image Build Targets"
	@echo "====================================="
	@echo ""
	@echo "Build all images:"
	@echo "  make docker-all                    - Build all Docker images"
	@echo ""
	@echo "Build individual images:"
	@echo "  make docker-protein-properties  - Protein properties"
	@echo "  make docker-zoogle             - Zoogle"
	@echo "  make docker-asteroid               - Asteroid"
	@echo "  make docker-bioservices            - Bioservices"
	@echo "  make docker-cialign                - CIAlign"
	@echo "  make docker-clipkit                - ClipKIT"
	@echo "  make docker-cogeqc                 - CoGeQC"
	@echo "  make docker-famsa                  - FAMSA"
	@echo "  make docker-fasttree               - FastTree"
	@echo "  make docker-generax                - GeneRax"
	@echo "  make docker-iqtree                 - IQ-TREE"
	@echo "  make docker-orthofinder            - OrthoFinder"
	@echo "  make docker-phylo-profiles         - Phylogenetic profiles"
	@echo "  make docker-rbase                  - R base image"
	@echo "  make docker-select-inflation       - MCL inflation selection"
	@echo "  make docker-witch                  - WITCH"
	@echo "  make docker-preprocess-proteomes   - Proteome preprocessing"
	@echo "  make docker-prequal                - PREQUAL"
	@echo ""
	@echo "Push images:"
	@echo "  make push-all                      - Push all images to Docker Hub"
	@echo "  make push-<image-name>             - Push specific image"
	@echo ""
	@echo "Other:"
	@echo "  make clean                         - Remove dangling Docker images"
	@echo ""
	@echo "Platform: $(DOCKER_PLATFORM)"
	@echo "Organization: $(DOCKER_ORG)"

# =============================================================================
# Build all images
# =============================================================================

docker-all: docker-protein-properties docker-zoogle \
	docker-asteroid docker-bioservices docker-cialign docker-clipkit \
	docker-cogeqc docker-famsa docker-fasttree docker-generax \
	docker-iqtree docker-orthofinder docker-phylo-profiles \
	docker-rbase docker-select-inflation docker-witch \
	docker-preprocess-proteomes docker-prequal
	@echo ""
	@echo "All Docker images built successfully!"

# =============================================================================
# Build individual images
# =============================================================================

# Images requiring root context
docker-protein-properties:
	@echo "Building $(PROTEIN_PROPERTIES_IMAGE):$(PROTEIN_PROPERTIES_TAG) from root context..."
	docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(PROTEIN_PROPERTIES_IMAGE):$(PROTEIN_PROPERTIES_TAG) \
		-f docker/protein_properties/Dockerfile \
		.
	@echo "Built successfully!"

docker-zoogle:
	@echo "Building $(ZOOGLE_IMAGE):$(ZOOGLE_TAG) from root context..."
	docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(ZOOGLE_IMAGE):$(ZOOGLE_TAG) \
		-f docker/zoogle/Dockerfile \
		.
	@echo "Built successfully!"

# Images built from their docker/ subdirectory
docker-asteroid:
	@echo "Building $(ASTEROID_IMAGE):$(ASTEROID_TAG)..."
	cd docker/asteroid && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(ASTEROID_IMAGE):$(ASTEROID_TAG) \
		.
	@echo "Built successfully!"

docker-bioservices:
	@echo "Building $(BIOSERVICES_IMAGE):$(BIOSERVICES_TAG)..."
	cd docker/bioservices && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(BIOSERVICES_IMAGE):$(BIOSERVICES_TAG) \
		.
	@echo "Built successfully!"

docker-cialign:
	@echo "Building $(CIALIGN_IMAGE):$(CIALIGN_TAG)..."
	cd docker/cialign && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(CIALIGN_IMAGE):$(CIALIGN_TAG) \
		.
	@echo "Built successfully!"

docker-clipkit:
	@echo "Building $(CLIPKIT_IMAGE):$(CLIPKIT_TAG)..."
	cd docker/clipkit && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(CLIPKIT_IMAGE):$(CLIPKIT_TAG) \
		.
	@echo "Built successfully!"

docker-cogeqc:
	@echo "Building $(COGEQC_IMAGE):$(COGEQC_TAG)..."
	cd docker/cogeqc && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(COGEQC_IMAGE):$(COGEQC_TAG) \
		.
	@echo "Built successfully!"

docker-famsa:
	@echo "Building $(FAMSA_IMAGE):$(FAMSA_TAG)..."
	cd docker/famsa && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(FAMSA_IMAGE):$(FAMSA_TAG) \
		.
	@echo "Built successfully!"

docker-fasttree:
	@echo "Building $(FASTTREE_IMAGE):$(FASTTREE_TAG)..."
	cd docker/fasttree && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(FASTTREE_IMAGE):$(FASTTREE_TAG) \
		.
	@echo "Built successfully!"

docker-generax:
	@echo "Building $(GENERAX_IMAGE):$(GENERAX_TAG)..."
	cd docker/generax && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(GENERAX_IMAGE):$(GENERAX_TAG) \
		.
	@echo "Built successfully!"

docker-iqtree:
	@echo "Building $(IQTREE_IMAGE):$(IQTREE_TAG)..."
	cd docker/iqtree && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(IQTREE_IMAGE):$(IQTREE_TAG) \
		.
	@echo "Built successfully!"

docker-orthofinder:
	@echo "Building $(ORTHOFINDER_IMAGE):$(ORTHOFINDER_TAG)..."
	cd docker/orthofinder && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(ORTHOFINDER_IMAGE):$(ORTHOFINDER_TAG) \
		.
	@echo "Built successfully!"

docker-phylo-profiles:
	@echo "Building $(PHYLO_PROFILES_IMAGE):$(PHYLO_PROFILES_TAG)..."
	cd docker/phylo_profiles && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(PHYLO_PROFILES_IMAGE):$(PHYLO_PROFILES_TAG) \
		.
	@echo "Built successfully!"

docker-rbase:
	@echo "Building $(RBASE_IMAGE):$(RBASE_TAG)..."
	cd docker/rbase && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(RBASE_IMAGE):$(RBASE_TAG) \
		.
	@echo "Built successfully!"

docker-select-inflation:
	@echo "Building $(SELECT_INFLATION_IMAGE):$(SELECT_INFLATION_TAG)..."
	cd docker/select_mcl_inflation_params && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(SELECT_INFLATION_IMAGE):$(SELECT_INFLATION_TAG) \
		.
	@echo "Built successfully!"

docker-witch:
	@echo "Building $(WITCH_IMAGE):$(WITCH_TAG)..."
	cd docker/witch && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(WITCH_IMAGE):$(WITCH_TAG) \
		.
	@echo "Built successfully!"

docker-preprocess-proteomes:
	@echo "Building $(PREPROCESS_PROTEOMES_IMAGE):$(PREPROCESS_PROTEOMES_TAG)..."
	cd docker/preprocess_proteomes && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(PREPROCESS_PROTEOMES_IMAGE):$(PREPROCESS_PROTEOMES_TAG) \
		.
	@echo "Built successfully!"

docker-prequal:
	@echo "Building $(PREQUAL_IMAGE):$(PREQUAL_TAG)..."
	cd docker/prequal && docker build \
		--platform $(DOCKER_PLATFORM) \
		-t $(PREQUAL_IMAGE):$(PREQUAL_TAG) \
		.
	@echo "Built successfully!"

# =============================================================================
# Push all images
# =============================================================================

push-all: push-protein-properties push-zoogle \
	push-asteroid push-bioservices push-cialign push-clipkit \
	push-cogeqc push-famsa push-fasttree push-generax \
	push-iqtree push-orthofinder push-phylo-profiles \
	push-rbase push-select-inflation push-witch \
	push-preprocess-proteomes push-prequal
	@echo ""
	@echo "All images pushed to Docker Hub successfully!"

# =============================================================================
# Push individual images
# =============================================================================

push-protein-properties: docker-protein-properties
	@echo "Pushing $(PROTEIN_PROPERTIES_IMAGE):$(PROTEIN_PROPERTIES_TAG)..."
	docker push $(PROTEIN_PROPERTIES_IMAGE):$(PROTEIN_PROPERTIES_TAG)
	@echo "Pushed successfully!"

push-zoogle: docker-zoogle
	@echo "Pushing $(ZOOGLE_IMAGE):$(ZOOGLE_TAG)..."
	docker push $(ZOOGLE_IMAGE):$(ZOOGLE_TAG)
	@echo "Pushed successfully!"

push-asteroid: docker-asteroid
	@echo "Pushing $(ASTEROID_IMAGE):$(ASTEROID_TAG)..."
	docker push $(ASTEROID_IMAGE):$(ASTEROID_TAG)
	@echo "Pushed successfully!"

push-bioservices: docker-bioservices
	@echo "Pushing $(BIOSERVICES_IMAGE):$(BIOSERVICES_TAG)..."
	docker push $(BIOSERVICES_IMAGE):$(BIOSERVICES_TAG)
	@echo "Pushed successfully!"

push-cialign: docker-cialign
	@echo "Pushing $(CIALIGN_IMAGE):$(CIALIGN_TAG)..."
	docker push $(CIALIGN_IMAGE):$(CIALIGN_TAG)
	@echo "Pushed successfully!"

push-clipkit: docker-clipkit
	@echo "Pushing $(CLIPKIT_IMAGE):$(CLIPKIT_TAG)..."
	docker push $(CLIPKIT_IMAGE):$(CLIPKIT_TAG)
	@echo "Pushed successfully!"

push-cogeqc: docker-cogeqc
	@echo "Pushing $(COGEQC_IMAGE):$(COGEQC_TAG)..."
	docker push $(COGEQC_IMAGE):$(COGEQC_TAG)
	@echo "Pushed successfully!"

push-famsa: docker-famsa
	@echo "Pushing $(FAMSA_IMAGE):$(FAMSA_TAG)..."
	docker push $(FAMSA_IMAGE):$(FAMSA_TAG)
	@echo "Pushed successfully!"

push-fasttree: docker-fasttree
	@echo "Pushing $(FASTTREE_IMAGE):$(FASTTREE_TAG)..."
	docker push $(FASTTREE_IMAGE):$(FASTTREE_TAG)
	@echo "Pushed successfully!"

push-generax: docker-generax
	@echo "Pushing $(GENERAX_IMAGE):$(GENERAX_TAG)..."
	docker push $(GENERAX_IMAGE):$(GENERAX_TAG)
	@echo "Pushed successfully!"

push-iqtree: docker-iqtree
	@echo "Pushing $(IQTREE_IMAGE):$(IQTREE_TAG)..."
	docker push $(IQTREE_IMAGE):$(IQTREE_TAG)
	@echo "Pushed successfully!"

push-orthofinder: docker-orthofinder
	@echo "Pushing $(ORTHOFINDER_IMAGE):$(ORTHOFINDER_TAG)..."
	docker push $(ORTHOFINDER_IMAGE):$(ORTHOFINDER_TAG)
	@echo "Pushed successfully!"

push-phylo-profiles: docker-phylo-profiles
	@echo "Pushing $(PHYLO_PROFILES_IMAGE):$(PHYLO_PROFILES_TAG)..."
	docker push $(PHYLO_PROFILES_IMAGE):$(PHYLO_PROFILES_TAG)
	@echo "Pushed successfully!"

push-rbase: docker-rbase
	@echo "Pushing $(RBASE_IMAGE):$(RBASE_TAG)..."
	docker push $(RBASE_IMAGE):$(RBASE_TAG)
	@echo "Pushed successfully!"

push-select-inflation: docker-select-inflation
	@echo "Pushing $(SELECT_INFLATION_IMAGE):$(SELECT_INFLATION_TAG)..."
	docker push $(SELECT_INFLATION_IMAGE):$(SELECT_INFLATION_TAG)
	@echo "Pushed successfully!"

push-witch: docker-witch
	@echo "Pushing $(WITCH_IMAGE):$(WITCH_TAG)..."
	docker push $(WITCH_IMAGE):$(WITCH_TAG)
	@echo "Pushed successfully!"

push-preprocess-proteomes: docker-preprocess-proteomes
	@echo "Pushing $(PREPROCESS_PROTEOMES_IMAGE):$(PREPROCESS_PROTEOMES_TAG)..."
	docker push $(PREPROCESS_PROTEOMES_IMAGE):$(PREPROCESS_PROTEOMES_TAG)
	@echo "Pushed successfully!"

push-prequal: docker-prequal
	@echo "Pushing $(PREQUAL_IMAGE):$(PREQUAL_TAG)..."
	docker push $(PREQUAL_IMAGE):$(PREQUAL_TAG)
	@echo "Pushed successfully!"

# =============================================================================
# Clean up
# =============================================================================

clean:
	@echo "Removing dangling Docker images..."
	docker image prune -f
	@echo "Cleanup complete!"
