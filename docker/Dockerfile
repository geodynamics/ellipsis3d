FROM ubuntu:jammy

RUN apt-get update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt-get install --yes \
    build-essential \
    git \
    cmake \
    automake \
    autoconf

RUN useradd \
  --create-home \
  ellipsis3d_user 

USER ellipsis3d_user

WORKDIR /home/ellipsis3d_user

RUN git clone 'https://github.com/geodynamics/ellipsis3d'; cd ellipsis3d; mkdir aux-config; autoreconf || true; automake --add-missing; ./configure; make;

ENV PATH="/home/ellipsis3d_user/ellipsis3d:${PATH}"
