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
# Conda Julia doesn't install the SSL certificates correctly
# RUN ln -s /etc/ssl/certs/ca-certificates.crt /opt/conda/share/julia/cert.pem

# Some dodgy dynamic linking conda doesn't set LD_LIBRARY_PATH
# ENV LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"

# Let's abandon Conda julia and install that separately
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_BASE_VERSION}/${JULIA_TARBALL}
RUN tar -xf $JULIA_TARBALL

# Enable our conda environment and pip-install pyjulia so we don't end up with two versions of Julia
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install julia
RUN python -c "import julia; julia.install(julia='${JULIA}')"

RUN $JULIA -e 'using Pkg; pkg"registry add https://github.com/JuliaRegistries/General"; pkg"registry add https://github.com/ACEsuit/ACEregistry"'

RUN PYTHON=$(which python) $JULIA -e 'using Pkg; Pkg.add(["ACE1", "ACE1x", "ASE", "JuLIP"])'

# Last working version:
RUN pip install git+https://gitlab.com/ase/ase.git@2481069f
# RUN pip install git+https://gitlab.com/ase/ase

RUN pip install git+https://github.com/casv2/pyjulip
RUN pip install git+https://github.com/ACEsuit/ACEHAL

# Sanity check that things are loading
# RUN PATH=$(dirname ${JULIA}):$PATH python -c "import pyjulip"

COPY --chown=$MAMBA_USER:$MAMBA_USER entrypoint.sh /tmp/entrypoint.sh
RUN chmod +x /tmp/entrypoint.sh

COPY --chown=$MAMBA_USER:$MAMBA_USER emt_training.py /tmp/emt_training.py
COPY --chown=$MAMBA_USER:$MAMBA_USER run_acehal.py /tmp/run_acehal.py

RUN echo "export PATH=$(dirname ${JULIA}):$PATH" >> /home/mambauser/.bashrc

ENTRYPOINT ["/tmp/entrypoint.sh"]
