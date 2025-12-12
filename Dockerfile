FROM julia:1.11.7
LABEL maintainer="Maarten Pronk <maarten.pronk@deltares.nl>"

RUN apt-get update && apt-get install -y \
    g++ git gcc && \
    rm -rf /var/lib/apt/lists/*
ADD . /app
WORKDIR /app/build/create_binaries/
RUN julia --project=/app/Wflow -e "using Pkg; Pkg.instantiate();" && \
    julia --project=/app/Wflow /app/utils/download_test_data.jl && \
    julia --project -e "using Pkg; Pkg.instantiate()" && \
    julia --project create_app.jl && \
    rm -rf /app/test/data/*

ENTRYPOINT [ "/app/build/create_binaries/wflow_bundle/bin/wflow_cli" ]
