FROM ubuntu:18.04

# Install prerequisites
RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" TZ="Europe/Berlin" \
    apt-get install -y \
    git \
    wget unzip \
    build-essential petsc-dev libeigen3-dev libxml2-dev libboost-all-dev \
    cmake \
    libffi-dev \
    python3 \
    vim

# make python point to python3 (needed for petsc build)
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 100

WORKDIR /workspace

# Clone opendihu
RUN git clone --branch docker-ubuntu16-18 --depth=1 https://github.com/maierbn/opendihu.git  
WORKDIR opendihu

# Build opendihu
RUN make release_without_tests; echo "done"; cat config.log

# Add aliases
ADD .bash_aliases /workspace
RUN echo ". /workspace/.bash_aliases " >> ~/.bashrc 

# Download input files
RUN cd examples/electrophysiology && wget https://zenodo.org/records/4705982/files/input.tgz && tar xf input.tgz
