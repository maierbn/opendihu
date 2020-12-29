# use the prerequisites image
FROM maierbn/opendihu-prerequisites:latest

# Set the working directory to /opendihu
WORKDIR /workspace

# checkout and prepare for opendihu
RUN apt-get update && \
  git clone https://github.com/maierbn/opendihu.git && \
  cd opendihu && \
  git checkout develop

# Build opendihu
RUN cd opendihu && make debug_without_tests release_without_tests; echo "done"; cat config.log

# Run system tests when the container launches
#CMD ["make", "system_testing"]
CMD ["make", "release"]

