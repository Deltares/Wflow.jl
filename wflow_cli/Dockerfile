FROM julia:1.8.2
LABEL maintainer="Maarten Pronk <maarten.pronk@deltares.nl>"

RUN apt-get update && apt-get install -y \
    g++ gcc \
    && rm -rf /var/lib/apt/lists/*
ADD . /app
WORKDIR /app
RUN julia contrib/build_steps.jl

ENTRYPOINT [ "/app/create_app/wflow_bundle/bin/wflow_cli" ]
