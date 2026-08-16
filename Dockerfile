FROM julia:1.12.7
LABEL maintainer="Maarten Pronk <maarten.pronk@deltares.nl>"

RUN apt-get update && apt-get install -y \
    g++ git gcc && \
    rm -rf /var/lib/apt/lists/*
ADD . /app
WORKDIR /app/build/create_binaries/
RUN julia --project=/app/Wflow -e "using Pkg; Pkg.instantiate();" && \
    julia --project -e "using Pkg; Pkg.instantiate()" && \
    julia --project create_app.jl

ENTRYPOINT [ "/app/build/create_binaries/wflow_bundle/bin/wflow_cli" ]
