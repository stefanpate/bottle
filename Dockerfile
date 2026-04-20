# Build layer
FROM python:3.12-bookworm AS builder

RUN pip install poetry==2.1.3

ENV POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1 \
    POETRY_CACHE_DIR=/tmp/poetry_cache

WORKDIR /app

COPY pyproject.toml poetry.lock ./

# Poetry doesn't like there to not be a README
RUN touch README.md

RUN poetry install --no-root --no-directory --only path_viewer && rm -rf $POETRY_CACHE_DIR

# Runtime layer
FROM python:3.12-slim-bookworm

ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:$PATH" \
    CASP_STUDY_ROOT=/data/processed \
    KNOWN_BIOCHEM_ROOT=/data/known \
    FEEDBACK_ROOT=/feedback \
    PYTHONPATH=/app

RUN apt-get update && apt-get install -y --no-install-recommends \
    libxrender1 \
    libxext6 \
    libexpat1 \
    rclone \
    sqlite3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=builder ${VIRTUAL_ENV} ${VIRTUAL_ENV}

# Code
COPY path_viewer ./path_viewer
COPY src ./src
COPY scripts ./scripts

# CASP and known data
COPY --from=casp . ${CASP_STUDY_ROOT}
COPY --from=known . ${KNOWN_BIOCHEM_ROOT}
COPY artifacts/plus.svg ./artifacts/plus.svg
COPY artifacts/arrow.svg ./artifacts/arrow.svg

RUN mkdir -p ${FEEDBACK_ROOT}

EXPOSE 8501

CMD ["streamlit", "run", "path_viewer/app.py", "--server.address=0.0.0.0", "--server.port=8501"]

