# How to use Docker for OpenDiHu

- Prerequisit: [install Docker](https://docs.docker.com/engine/install/ubuntu/)
- Build a provided OpenDiHu container:
    1. Move to the directory where `Dockerfile` is located
    2. Execute `docker build -t workspace .`
        > [!NOTE]  
        > - You can add the flag `--progress=plain` for more verbosity
        > - You can add the flag `--no-cache` for a clean build
- Run the container you just built: `docker run -i -t workspace`
