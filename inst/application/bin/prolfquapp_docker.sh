#!/bin/bash
# This script is used to run a single command in the prolfquapp docker image.
# - The docker image is pulled if it does not exist locally.
# - The current directory will be mounted to /work and set as the current working directory for the command being executed.
# Limitations
# - Currently, specified paths should be relative to the current working directory and not go outside of it.

set -euo pipefail

# Default values for arguments
IMAGE_VERSION="0.0.5"
IMAGE_REPO="docker.io/prolfqua/prolfquapp"

# Function to print usage/help message
usage() {
  echo "Usage: $0 [--image-version VERSION] [--image-repo REPO] [container_args]"
  echo "Options:"
  echo "  --image-version   Specify the version of the image to use (default: $IMAGE_VERSION)"
  echo "  --image-repo      Specify the repository of the image (default: $IMAGE_REPO)"
  exit 1
}

if [ "$#" -eq 0 ]; then
  usage
fi

# Argument parsing
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --image-version)
            IMAGE_VERSION="$2"
            shift 2
            ;;
        --image-repo)
            IMAGE_REPO="$2"
            shift 2
            ;;
        --help)
            usage
            ;;
        --) # Stop parsing options
            shift
            break
            ;;
        -*)
            echo "Unknown option: $1"
            usage
            ;;
        *)
            break
            ;;
    esac
done

# Remaining arguments after known options are treated as container_args
CONTAINER_ARGS=("$@")

# Function to run the container
run() {
    local image_version="$1"
    local image_repo="$2"
    shift 2
    local container_args=("$@")

    # Check if podman is installed, otherwise use docker
    if command -v podman > /dev/null; then
        DOCKER="podman"
    else
        DOCKER="docker"
    fi

    local image="${image_repo}:${image_version}"
    echo "Using $DOCKER to run $image"

    # Check if the image is available locally
    if ! $DOCKER image inspect "$image" > /dev/null 2>&1; then
        echo "Image $image not found locally. Pulling..."
        $DOCKER pull "$image"
    else
        echo "Image $image already exists locally."
        echo "If you want to update the image to the latest version, run the following command:"
        echo "$DOCKER pull $image"
    fi

    # Run the container with current user UID/GID and mount the current working directory
    $DOCKER run --user "$(id -u):$(id -g)" --rm -it --mount "type=bind,source=$(pwd),target=/work" -w /work "$image" "${container_args[@]}"
}

# Call the run function with the parsed arguments
run "$IMAGE_VERSION" "$IMAGE_REPO" "${CONTAINER_ARGS[@]}"
