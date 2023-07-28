FROM julia:1.9.2
LABEL maintainer="Maarten Pronk <maarten.pronk@deltares.nl>"

RUN apt-get update && apt-get install -y \
    g++ gcc \
    && rm -rf /var/lib/apt/lists/*
ADD . /app
WORKDIR /app/build/create_binaries/
RUN julia --project -e "using Pkg; Pkg.instantiate()"
RUN julia --project create_app.jl

ENTRYPOINT [ "/app/build/create_binaries/wflow_bundle/bin/wflow_cli" ]
