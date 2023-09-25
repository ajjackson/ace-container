FROM mambaorg/micromamba
COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# Conda Julia doesn't install the SSL certificates correctly
RUN ln -s /etc/ssl/certs/ca-certificates.crt /opt/conda/share/julia/cert.pem

# Some dodgy dynamic linking conda doesn't set LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"

RUN python -c "import julia; julia.install()"

RUN julia -e 'using Pkg; pkg"registry add https://github.com/JuliaRegistries/General"; pkg"registry add https://github.com/ACEsuit/ACEregistry"'
RUN julia -e 'using Pkg; Pkg.add(["ACE1", "ACE1x", "ASE", "JuLIP"])'

RUN pip install git+https://github.com/casv2/pyjulip
RUN pip install git+https://gitlab.com/ase/ase
RUN pip install git+https://github.com/ACEsuit/ACEHAL
RUN python -c "import pyjulip"
