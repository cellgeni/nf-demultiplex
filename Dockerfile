FROM python:3.10

# use bash shell
SHELL ["/bin/bash", "-c"]

# non-interactive mode
# make apt-get have zero interaction while installing or upgrading
ENV DEBIAN_FRONTEND=noninteractive

# trying to match versions of tools and packages to official image

# OS dependencies
RUN apt-get update -y && \
    apt-get install -y wget build-essential \
      libncurses5-dev zlib1g-dev \
      libbz2-dev liblzma-dev procps

# minimap2
ARG MINIMAP2="2.26"
RUN wget -q -O minimap2.tar.bz2 https://github.com/lh3/minimap2/releases/download/v${MINIMAP2}/minimap2-${MINIMAP2}_x64-linux.tar.bz2 && \
    mkdir -p /opt/minimap2 && \
    tar xjvf minimap2.tar.bz2 -C /opt/minimap2 --strip-components 1
ENV PATH="/opt/minimap2:$PATH"

# bedtools  
ARG BEDTOOLS2="2.28.0"
RUN wget -q -O bedtools2.tar.gz https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS2}/bedtools-${BEDTOOLS2}.tar.gz && \
    mkdir -p /opt/bedtools2 && \
    tar zxvf bedtools2.tar.gz  -C /opt/bedtools2 --strip-components 1 && \
    cd /opt/bedtools2 && make
ENV PATH="/opt/bedtools2/bin:$PATH"

# Rust
ENV CARGO_HOME=/opt/rust/.cargo
ENV RUSTUP_HOME=/opt/rust/.cargo 
RUN wget -qO- https://sh.rustup.rs | bash -s -- -y && \
    . "/opt/rust/.cargo/env" && \
    rustup default stable
ENV PATH="/opt/rust/.cargo/bin:$PATH"

# souporcell
ARG SOUPORCELL="2.5"
RUN wget -q -O "souporcell.tar.gz" "https://github.com/wheaton5/souporcell/archive/refs/tags/v${SOUPORCELL}.tar.gz" && \
    mkdir -p /opt/souporcell && \
    tar zxvf souporcell.tar.gz  -C /opt/souporcell --strip-components 1 && \
    cd /opt/souporcell/troublet && \
    cargo build --release && \
    cd /opt/souporcell/souporcell && \
    cargo build --release
ENV PATH="/opt/souporcell:/opt/souporcell/troublet/target/release:$PATH"

# python dependencies (trying to use exactly the same as the ones in the official image)
RUN pip install --upgrade pip && \
    # to prevent 'error in PyVCF setup command: use_2to3 is invalid' need <58.0.0 
    pip install --no-cache-dir "setuptools==57.5.0" && \
    pip install --no-cache-dir Cython && \
    pip install --no-cache-dir pysam pyvcf "numpy==1.24.3" scipy pyfaidx && \
    pip install --no-cache-dir "pystan==2.19.1.1"
      
#pysam pyvcf "numpy<1.23.0" scipy "Cython>=0.22" pyfaidx && \

# HTSlib/samtools/bcftools
ARG SAMTOOLS="1.9"
RUN wget -q -O htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/${SAMTOOLS}/htslib-${SAMTOOLS}.tar.bz2 && \
    mkdir -p /opt/histlib && \
    tar xvfj htslib.tar.bz2 -C /opt/histlib --strip-components 1 && \
    cd /opt/histlib && ./configure && make && make install
RUN wget -q -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLS}/samtools-${SAMTOOLS}.tar.bz2 && \
    mkdir -p /opt/samtools && \
    tar xvfj samtools.tar.bz2 -C /opt/samtools --strip-components 1 && \
    cd /opt/samtools && ./configure && make && make install
RUN wget -q -O bcftools.tar.bz2 https://github.com/samtools/bcftools/releases/download/${SAMTOOLS}/bcftools-${SAMTOOLS}.tar.bz2 && \
    mkdir -p /opt/bcftools && \
    tar xvfj bcftools.tar.bz2 -C /opt/bcftools --strip-components 1 && \
    cd /opt/bcftools && ./configure && make && make install

# freebayes
ARG FREEBAYES="1.3.6"
# This package doesnt keep consistent naming so when updating may need to change
RUN wget -q -O /opt/freebayes https://github.com/freebayes/freebayes/releases/download/v${FREEBAYES}/freebayes-${FREEBAYES}-linux-amd64-static.gz && \
    chmod +x /opt/freebayes

# Vartrix
ARG VARTRIX="1.1.22"
RUN wget -q -O /opt/vartrix https://github.com/10XGenomics/vartrix/releases/download/v${VARTRIX}/vartrix_linux && \
    chmod +x /opt/vartrix
ENV PATH="/opt:$PATH"

# remove the entrypoint from the parent container
ENTRYPOINT ["bash"]
