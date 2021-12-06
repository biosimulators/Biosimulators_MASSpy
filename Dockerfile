# Base OS
FROM python:3.9-slim-buster

ARG VERSION="0.0.5"
ARG SIMULATOR_VERSION="0.1.5"

# metadata
LABEL \
    org.opencontainers.image.title="MASSpy" \
    org.opencontainers.image.version="${SIMULATOR_VERSION}" \
    org.opencontainers.image.description="Tool for kinetic simulation of metabolic networks." \
    org.opencontainers.image.url="https://masspy.readthedocs.io/" \
    org.opencontainers.image.documentation="https://masspy.readthedocs.io/" \
    org.opencontainers.image.source="https://github.com/SBRG/MASSpy" \
    org.opencontainers.image.authors="BioSimulators Team <info@biosimulators.org>" \
    org.opencontainers.image.vendor="BioSimulators Team" \
    org.opencontainers.image.licenses="SPDX:MIT" \
    \
    base_image="python:3.9-slim-buster" \
    version="${VERSION}" \
    software="MASSpy" \
    software.version="${SIMULATOR_VERSION}" \
    about.summary="Tool for kinetic simulation of metabolic networks." \
    about.home="https://masspy.readthedocs.io/" \
    about.documentation="https://masspy.readthedocs.io/" \
    about.license_file="https://github.com/SBRG/MASSpy/blob/master/LICENSE" \
    about.license="SPDX:MIT" \
    about.tags="kinetic modeling,dynamical simulation,systems biology,biochemical networks,MASSpy,SED-ML,COMBINE,OMEX,BioSimulators" \
    maintainer="BioSimulators Team <info@biosimulators.org>"

# Install requirements for tellurium
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        libxml2 \
        libncurses5 \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Create COBRApy cache
RUN mkdir -p /.cache/cobrapy

# Copy code for command-line interface into image and install it
COPY . /root/Biosimulators_MASSpy
RUN pip install /root/Biosimulators_MASSpy \
    && rm -rf /root/Biosimulators_MASSpy
ENV VERBOSE=0 \
    MPLBACKEND=PDF
RUN mkdir -p /.config/matplotlib \
    && mkdir -p /.cache/matplotlib \
    && chmod ugo+rw /.config/matplotlib \
    && chmod ugo+rw /.cache/matplotlib

# Entrypoint
ENTRYPOINT ["biosimulators-masspy"]
CMD []
