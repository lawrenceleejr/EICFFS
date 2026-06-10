FROM mambaorg/micromamba:1.5.3

LABEL maintainer="EIC FFS Study"
LABEL description="Docker container for EIC FFS (Frame-dependent Fragmentation Shift) phenomenology study"

USER root

# Install system dependencies for build tools
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        libgsl-dev \
        wget \
        git \
    && rm -rf /var/lib/apt/lists/*

# Create working directory
RUN mkdir -p /work && chown $MAMBA_USER:$MAMBA_USER /work

USER $MAMBA_USER

# Copy and install conda environment
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

# Set PATH so conda environment is active
ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /work

# Copy all source files
COPY --chown=$MAMBA_USER:$MAMBA_USER generate_events.py /work/
COPY --chown=$MAMBA_USER:$MAMBA_USER analyze_jets.py /work/
COPY --chown=$MAMBA_USER:$MAMBA_USER make_results.py /work/
COPY --chown=$MAMBA_USER:$MAMBA_USER make_paper_plots.py /work/
COPY --chown=$MAMBA_USER:$MAMBA_USER convert_hepmc.py /work/
COPY --chown=$MAMBA_USER:$MAMBA_USER run_herwig.sh /work/
COPY --chown=$MAMBA_USER:$MAMBA_USER herwig/ /work/herwig/
COPY --chown=$MAMBA_USER:$MAMBA_USER analyze_events.py /work/
COPY --chown=$MAMBA_USER:$MAMBA_USER make_plots.py /work/
COPY --chown=$MAMBA_USER:$MAMBA_USER utils/ /work/utils/
COPY --chown=$MAMBA_USER:$MAMBA_USER run.sh /work/

RUN chmod +x /work/run.sh

# Default: run the full pipeline
CMD ["/work/run.sh"]
