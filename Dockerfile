FROM mambaorg/micromamba
COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# Conda Julia doesn't install the SSL certificates correctly
RUN ln -s /etc/ssl/certs/ca-certificates.crt /opt/conda/share/julia/cert.pem

# Help Julia find Conda
ARG CONDA_JL_HOME=/opt/conda/envs/base
ARG PYTHON=/opt/conda/bin/python3
ARG PYCALL_JL_RUNTIME_PYTHON=/opt/conda/bin/python3
ARG CONDA_JL_CONDA_EXE=/opt/conda/bin/conda
ENV LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"

RUN julia -e 'using Pkg; pkg"registry add https://github.com/JuliaRegistries/General"; pkg"registry add https://github.com/ACEsuit/ACEregistry"'

# RUN julia -e 'using Pkg; Pkg.add("Conda")'
# RUN julia -e 'using Pkg; Pkg.add("PyCall"); Pkg.build("PyCall")'
# RUN julia -e 'using Pkg; pkg"add ACE1, ACE1x, ASE, JuLIP"'

RUN python -m pip install julia==0.6.1

# RUN pip install git+https://github.com/ACEsuit/ACEHAL

RUN python -c "import julia; julia.install()"

RUN julia -e 'using Pkg; pkg"add JuLIP"'
RUN pip install git+https://github.com/casv2/pyjulip
RUN pip install git+https://gitlab.com/ase/ase
RUN python -c "import pyjulip"
