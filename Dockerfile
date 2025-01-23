#
# --- python build layer
#
FROM synbiorox/base:0.1.0 AS builder

# TODO look into removing this dep
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libopenblas-dev \
    libxrender-dev \
    libxext-dev

ARG POETRY_VERSION=1.8

RUN pip install poetry==$POETRY_VERSION \
    && poetry self add poetry-plugin-bundle;

ENV VIRTUAL_ENV=/app/.venv \
    POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1 \
    POETRY_CACHE_DIR=/tmp/poetry_cache

WORKDIR /app

COPY src ./src
COPY README.md pyproject.toml poetry.lock ./
# RUN touch README.md

RUN poetry bundle venv $VIRTUAL_ENV \
    --without dev \
    --verbose \
    --no-cache \
    && rm -rf $POETRY_CACHE_DIR;

#
# --- runtime layer
#
FROM synbiorox/base:0.1.0 AS runtime


# FIXME using pixi with conda bin deps would simplify
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libopenblas-dev \
    libxrender-dev \
    libxext-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*;

ARG NB_USER=bottler
ARG NB_UID=1000
ARG NB_GROUP=b1039
ARG NB_GID=111111039

ENV VIRTUAL_ENV=/app/.venv \
    JAVA_HOME=/opt/java/openjdk \
    PATH="/app/.venv/bin:$PATH" \
    USER=${NB_USER} \
    NB_UID=${NB_UID} \
    HOME=/home/${NB_USER}

ENV BOTTLE_EXPANSION_ASSETS=${HOME}/assets \
    PATH="${JAVA_HOME}/bin:${PATH}"

RUN groupadd -g ${NB_GID} ${NB_GROUP} \
    && adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    --gid ${NB_GID} \
    ${NB_USER}

# copy python virualenv from the build layer
COPY --from=builder ${VIRTUAL_ENV} ${VIRTUAL_ENV}

# FIXME this is a rather clunky way of doing this
# copy project assets (provided from host)
WORKDIR ${BOTTLE_EXPANSION_ASSETS}
COPY --from=assets ./found_paths.json .
COPY --from=assets ./known_reactions.json .
COPY --from=assets ./predicted_reactions.json .
COPY --from=assets ./svgs ./svgs

# copy project source code
WORKDIR ${HOME}
COPY notebooks/ ${HOME}/notebooks
COPY voila.json ${HOME}/voila.json

# change ownership of project files to root
RUN chown -R ${NB_UID}:${NB_GID} ${HOME} \
    && chmod -R go+rX $HOME;

# set the user the default for runtime
USER ${NB_USER}

ENV IPYNB_PATH_VIEWER=notebooks/path_viewer.ipynb
RUN jupyter trust ${IPYNB_PATH_VIEWER}

# run the voila server with the TOKEN set and open to path_viewer
CMD ["sh", "-c", "voila $IPYNB_PATH_VIEWER --port=8888 --Voila.ip=0.0.0.0 --token --no-browser"]

# Jupyter testing
#ENV JUPYTER_TOKEN=dolalay
#CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser"]
