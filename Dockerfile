# --- build layer
FROM python:3.12-bookworm AS builder

# TODO look into removing this dep
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libopenblas-dev \
    libxrender-dev \
    libxext-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ARG POETRY_VERSION=2.1.0

RUN pip install poetry==$POETRY_VERSION \
    && poetry self add poetry-plugin-bundle

ENV VIRTUAL_ENV=/app/.venv \
    POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1 \
    POETRY_CACHE_DIR=/tmp/poetry_cache

WORKDIR /app

COPY src ./src
COPY README.md pyproject.toml poetry.lock ./
RUN touch README.md

RUN poetry bundle venv $VIRTUAL_ENV \
    --without dev \
    --verbose \
    --no-cache \
    && rm -rf $POETRY_CACHE_DIR

# --- runtime layer
FROM python:3.12-slim-bookworm AS runtime

# FIXME pixi using conda bin deps should fix
# potentially other ways of fixing this as well
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libopenblas-dev \
    libxrender-dev \
    libxext-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ARG NB_USER=jovyan
ARG NB_UID=1000

ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:$PATH" \
    USER=${NB_USER} \
    NB_UID=${NB_UID} \
    HOME=/home/${NB_USER}

ENV BOTTLE_EXPANSION_ASSETS=${HOME}/assets

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

COPY --from=builder ${VIRTUAL_ENV} ${VIRTUAL_ENV}


# copy project assets
WORKDIR ${BOTTLE_EXPANSION_ASSETS}
COPY --from=assets ./found_paths.parquet .
COPY --from=assets ./predicted_reactions.parquet .
COPY --from=assets ./svgs ./svgs

# copy project source code
WORKDIR ${HOME}
COPY notebooks/ ${HOME}/notebooks
COPY voila.json ${HOME}/voila.json
COPY ./artifacts/known/known_reactions.parquet ${HOME}/known/known_reactions.parquet
COPY ./artifacts/known/known_enzymes.parquet ${HOME}/known/known_enzymes.parquet

# change ownership of project files to root
RUN chown -R ${NB_UID} ${HOME}
# set the user the default for runtime
USER ${NB_USER}

ENV IPYNB_PATH_VIEWER=notebooks/path_viewer.ipynb
RUN jupyter trust ${IPYNB_PATH_VIEWER}

# run the voila server with the TOKEN set and open to path_viewer
CMD ["sh", "-c", "voila $IPYNB_PATH_VIEWER --port=8888 --Voila.ip=0.0.0.0 --token --no-browser"]

# Jupyter testing
#ENV JUPYTER_TOKEN=dolalay
#CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser"]