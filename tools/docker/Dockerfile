FROM ubuntu:22.04

# Install prerequisites
RUN apt update -y && apt install -y \
git wget vim \
build-essential cmake python3-dev python3-numpy python-is-python3 \
libeigen3-dev libxml2-dev libboost-all-dev petsc-dev 

WORKDIR /workspace

# Clone opendihu
RUN git clone --branch docker-with-precice-v2.5 --depth=1 https://github.com/maierbn/opendihu.git  
WORKDIR opendihu

# Build opendihu
RUN make release_without_tests; echo "done"; cat config.log

# Add aliases
ADD .bash_aliases /workspace
RUN echo ". /workspace/.bash_aliases " >> ~/.bashrc 

