#!/usr/bin/env zsh

# Project-specific variables
name="devcontainer"
version="latest"

# Stop on failure
set -e -o pipefail

# Set log file
exec > >(tee build.log) 2>&1

# Set platform
case $1 in
native)
    platform=$(uname -m)
    ;;
*)
    platform="amd64"
    ;;
esac

# Singularity image name
img="${name}_${version}_${platform}_$(date "+%Y_%m_%d_%H_%M_%S").sif"

# Build container with Docker
docker-compose build --progress=plain
#  --platform linux/${platform} -t ${name}:${version} .

# Convert Docker image to Singularity image
unset DOCKER_HOST
#export DOCKER_HOST=unix:///var/run/docker.sock
#export DOCKER_HOST=tcp://docker-desktop:2375
echo $DOCKER_HOST
docker run -v /var/run/docker.sock:/var/run/docker.sock\
 -v "$PWD:/output" --privileged -t --rm singularityware/docker2singularity\
 -n ${img} ${name}:${version}
mv ${img%.sif}.simg ${img}

# Change Singularity image permissions
if [[ $(uname -s) == "Linux" ]]; then
  sudo chown $USER:$USER $img
fi

# Update module file with new Singularity image name
sed -i'' -e "s/^local img_name.*/local img_name      = \'${img}\'/g"\
 ${name}.lua

