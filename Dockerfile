FROM mambaorg/micromamba
COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

USER root
RUN apt-get update && apt-get install -y wget git

ENV JULIA_BASE_VERSION=1.9
ENV JULIA_VERSION=1.9.1
ENV JULIA=/tmp/julia-${JULIA_VERSION}/bin/julia
ENV JULIA_TARBALL=julia-${JULIA_VERSION}-linux-x86_64.tar.gz

USER mambauser
RUN echo "export PATH=$(dirname ${JULIA}):\$PATH" >> /home/mambauser/.bashrc

# Let's abandon Conda julia and install from recommended package
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_BASE_VERSION}/${JULIA_TARBALL}
RUN tar -xf $JULIA_TARBALL

# Enable our conda environment and pip-install pyjulia so we don't end up with two versions of Julia
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install julia
RUN python -c "import julia; julia.install(julia='${JULIA}')"

RUN $JULIA -e 'using Pkg; pkg"registry add https://github.com/JuliaRegistries/General"; pkg"registry add https://github.com/ACEsuit/ACEregistry"'

RUN PYTHON=$(which python) $JULIA -e 'using Pkg; Pkg.add(["ACE1", "ACE1x", "ASE", "JuLIP"])'

# Last working version of ASE:
RUN pip install git+https://gitlab.com/ase/ase.git@2481069f
# RUN pip install git+https://gitlab.com/ase/ase

RUN pip install git+https://github.com/casv2/pyjulip
RUN pip install git+https://github.com/ACEsuit/ACEHAL

COPY --chown=$MAMBA_USER:$MAMBA_USER emt_training.py /tmp/emt_training.py
COPY --chown=$MAMBA_USER:$MAMBA_USER run_acehal.py /tmp/run_acehal.py
