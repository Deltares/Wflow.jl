FROM julia:1.9.2
LABEL maintainer="Maarten Pronk <maarten.pronk@deltares.nl>"

RUN apt-get update && apt-get install -y \
    g++ gcc \
    && rm -rf /var/lib/apt/lists/*
ADD . /app
WORKDIR /app
RUN julia build/create_binaries/create_app.jl

ENTRYPOINT [ "/app/create_app/wflow_bundle/bin/wflow_cli" ]
