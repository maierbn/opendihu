FROM ubuntu:22.04

# Install prerequisites
RUN apt-get update && \
    apt-get install -y \
    git \
    wget unzip \
    build-essential petsc-dev bison flex libeigen3-dev libxml2-dev libboost-all-dev \
    cmake \
    libffi-dev \
    vim

# git -> Dockerfile
# wget unzip -> Download dependencies
# build-essential petsc-dev bison flex libeigen3-dev libxml2-dev libboost-all-dev -> PETSc
# cmake -> base64
# libffi-dev -> Solves error at run time: ModuleNotFoundError: No module named '_ctypes’
# vim -> optional

# make `python` point to python3 (needed for some scripts that call python instead of python3)
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 100

WORKDIR /workspace

# Clone opendihu
RUN git clone --branch develop --depth=1 https://github.com/maierbn/opendihu.git  
WORKDIR opendihu

# Set environment variable such that OpenMPI works inside docker container (run as root)
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Build opendihu, debug and release target for core, build release unit tests and run tests, build all examples
RUN make; echo "done"; cat config.log

# Add aliases
ADD .bash_aliases /workspace
RUN echo ". /workspace/.bash_aliases " >> ~/.bashrc 

# Download input files
RUN cd examples/electrophysiology && wget https://zenodo.org/records/4705982/files/input.tgz && tar xf input.tgz