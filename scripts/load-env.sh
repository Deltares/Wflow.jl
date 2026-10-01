#!/usr/bin/env bash

if [[ -f "${PIXI_PROJECT_ROOT}/.env" ]]; then
    eval "$(dotenv -f "${PIXI_PROJECT_ROOT}/.env" list --format=export)"
fi
