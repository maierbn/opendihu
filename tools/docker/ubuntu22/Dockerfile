FROM ubuntu:22.04

# Install prerequisites
RUN apt-get update && apt-get install -y \
    git \
    wget unzip \
    build-essential petsc-dev libeigen3-dev libxml2-dev libboost-all-dev \
    cmake \
    libffi-dev \
    vim

# git -> Dockerfile
# wget unzip -> Download dependencies
# build-essential petsc-dev libeigen3-dev libxml2-dev libboost-all-dev -> PETSc
# cmake -> base64
# libffi-dev -> Solves error at run time: ModuleNotFoundError: No module named '_ctypes’
# vim -> optional

WORKDIR /workspace

# Clone opendihu
RUN git clone --branch docker-with-precice-v2.5 --depth=1 https://github.com/maierbn/opendihu.git  
WORKDIR opendihu

# Build opendihu
RUN make release_without_tests; echo "done"; cat config.log

# Add aliases
ADD .bash_aliases /workspace
RUN echo ". /workspace/.bash_aliases " >> ~/.bashrc 

# Download input files
RUN cd examples/electrophysiology && wget https://zenodo.org/records/4705982/files/input.tgz && tar xf input.tgz