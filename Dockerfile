FROM ubuntu:18.04

RUN mkdir -p /tmp/python-installation
RUN mkdir -p /software/msasa

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

RUN apt-get --yes -qq update \
    && apt-get --yes -qq upgrade \
    && apt-get --yes -qq install \
        bzip2 \
        cmake \
        cpio \
        curl \
        g++ \
        gcc \
        gfortran \
        git \
        gosu \
        libblas-dev \
        liblapack-dev \
        libopenmpi-dev \
        openmpi-bin \
        python3-dev \
        python3-pip \
        virtualenv \
        wget \
        zlib1g-dev \
        vim       \
        htop      \
        libncurses5-dev \
        libgdbm-dev \
        build-essential \
        libnss3-dev \
        libssl-dev \
        libreadline-dev \
        libffi-dev \
        libsqlite3-dev \
        libbz2-dev \
        patchelf \
        ccache
    #     # time \
    # && apt-get --yes -qq clean \
    # && rm -rf /var/lib/apt/lists/*

RUN /usr/sbin/update-ccache-symlinks
RUN export PATH="/usr/lib/ccache:$PATH"

WORKDIR /tmp/python-installation
RUN wget https://www.python.org/ftp/python/3.10.6/Python-3.10.6.tgz && tar -xf Python-3.10.6.tgz

WORKDIR /tmp/python-installation/Python-3.10.6
RUN ./configure --enable-optimizations
RUN make -j 8 && make altinstall

RUN python3.10 --version
# RUN python3.10 -m pip install --upgrade pip

WORKDIR /software/msasa
COPY ./requirements.txt /software/msasa/requirements.txt
RUN python3.10 -m pip install --no-input -r /software/msasa/requirements.txt

COPY ./src /software/msasa/src

RUN python3.10 -m pip install ordered-set nuitka
RUN python3.10 -m pip install pyqt5
RUN python3.10 -m nuitka \
    --static-libpython=yes \
    --follow-imports \
    --plugin-enable=numpy \
    --plugin-enable=pyqt5 \
    --assume-yes-for-downloads \
    --standalone \
    --lto=no \
    --show-progress \
    --show-memory \
    --jobs=12 \
    --show-modules src/msa.py

ENTRYPOINT [ "/software/msasa/msa.dist/msa" ]
# ENTRYPOINT [ "python3.10", "/software/msasa/src/msa.py" ]
