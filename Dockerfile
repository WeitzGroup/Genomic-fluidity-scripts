##########################################################
# Dockerfile to build a Genomic-fluidity-scripts container image
# Based on Ubuntu
############################################################

# Set the base image to the BioPerl prebuilt prerequisites image

FROM bioperl/bioperl-deps

# File Author / Maintainer
MAINTAINER Nelly Selem <nselem84@gmail.com>

# First Bioperl Modules copied from Hilmar Lapp <hlapp@drycafe.net> bioperl
# Install modules recommended by BioPerl.
# We can't include Bio::ASN1::EntrezGene here yet, because it declares
# a dependency on BioPerl, thus pulling in BioPerl first.
RUN cpanm \
  Bio::Phylo

# -------------------------------------------------------------
# Install BioPerl from GitHub current master branch.
#
# This is the actual installation step of BioPerl itself :-)
# -------------------------------------------------------------
RUN cpanm -v \
  https://github.com/bioperl/bioperl-live/archive/master.tar.gz

# Install modules recommended by BioPerl that depend on BioPerl.
# (Don't ask. See above.)
RUN cpanm \
  Bio::ASN1::EntrezGene


## Cloning fluidity scripts Weitz Group
RUN apt-get update && apt-get install -y git
RUN git clone https://github.com/WeitzGroup/Genomic-fluidity-scripts

## Bioperl module
RUN cpanm \
Bio::Tools::Run::StandAloneBlast

# Installing blast legacy
RUN apt-get install -y blast2

RUN mkdir /usr/src/fluidity
COPY . /usr/src/fluidity
ENV PATH $PATH:/Genomic-fluidity-scripts:
WORKDIR /Genomic-fluidity-scripts 
RUN chmod +x *
WORKDIR /usr/src/fluidity 

