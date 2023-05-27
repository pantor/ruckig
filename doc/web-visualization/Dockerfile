# https://hub.docker.com/r/emscripten/emsdk

# Point at any base image that you find suitable to extend.
FROM emscripten/emsdk:3.1.39

# Install required tools that are useful for your project i.e. ninja-build
RUN apt update && apt install -y ninja-build cmake bash
