# use ubuntu 20.04 as parent image
FROM ubuntu:20.04

# Set the working directory to /workspace
WORKDIR /workspace

# Copy the current directory contents into the container
ADD . /workspace

# Install prerequisites
RUN apt-get update && \
  DEBIAN_FRONTEND="noninteractive" TZ="Europe/Berlin" \
  apt-get install -y \
    git make cmake apt-utils software-properties-common libopenmpi-dev libx11-* zlib1g-dev libssl-dev libffi-dev bison flex \
    python-is-python3 python3-pip libboost-filesystem-dev libboost-log-dev libboost-program-options-dev libboost-system-dev libboost-thread-dev libboost-test-dev \
    vim tig valgrind gdb wget geany geany-plugins-* meld

# Load bash aliases, which is useful for interactive sessions
ADD .bash_aliases /workspace
RUN echo ". /workspace/.bash_aliases " >> ~/.bashrc && \
echo "set tabstop=2" >> ~/.vimrc &&  \
echo "set background=dark" > ~/.vimrc && \
echo "set syntax=python" > ~/.vimrc && \
echo "set ignorecase" > ~/.vimrc && \
echo "set smartcase" > ~/.vimrc && \
echo "set number" > ~/.vimrc && \
echo "set wrap!" > ~/.vimrc && \
echo "au BufRead,BufNewFile *.tpp set filetype=cpp" > ~/.vimrc && \
echo "set -g terminal-overrides 'xterm*:smcup@:rmcup@'" > ~/.tmux.conf && \
echo "set-option -g history-limit 3000000" >> ~/.tmux.conf

