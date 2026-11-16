# microbiomeutils container image.
#
# Multi-stage build: scikit-bio 0.5.9 (and its dependency hdmedians) have no
# prebuilt wheels for Python 3.11, so they are compiled in a builder stage and
# only the resulting wheels are installed in the slim runtime image.

ARG PYTHON_IMAGE=python:3.11-slim-bookworm

# --- Builder: compile wheels for all pinned dependencies -------------------
FROM ${PYTHON_IMAGE} AS builder

RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt constraints.txt /tmp/
# PIP_CONSTRAINT also applies to the isolated build environments, so the
# numpy/Cython used to compile scikit-bio and hdmedians are pinned too.
RUN PIP_CONSTRAINT=/tmp/constraints.txt \
    pip wheel --no-cache-dir -r /tmp/requirements.txt -w /wheels

# --- Runtime -------------------------------------------------------------
FROM ${PYTHON_IMAGE}

ARG VERSION=dev
LABEL org.opencontainers.image.title="microbiomeutils" \
      org.opencontainers.image.description="Beta diversity, PCoA, EMPeror plots and taxonomic summaries from a feature table" \
      org.opencontainers.image.source="https://github.com/jtremblay/microbiomeutils" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.licenses="MIT"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    MPLCONFIGDIR=/tmp/matplotlib

# procps provides ps, which workflow managers such as Nextflow need inside the
# container to collect task metrics (CPU, memory).
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt constraints.txt /tmp/
RUN --mount=type=bind,from=builder,source=/wheels,target=/wheels \
    pip install --no-cache-dir --no-index --find-links=/wheels \
        -r /tmp/requirements.txt -c /tmp/constraints.txt \
    && rm /tmp/requirements.txt /tmp/constraints.txt

# microbiomeutils.py imports nrc relative to its own directory: keep them side by side.
COPY microbiomeutils.py /opt/microbiomeutils/
COPY nrc/ /opt/microbiomeutils/nrc/
COPY tests/ /opt/microbiomeutils/tests/

RUN useradd --create-home --uid 1000 --user-group microbiome

USER microbiome
WORKDIR /data

ENTRYPOINT ["python3", "/opt/microbiomeutils/microbiomeutils.py"]
CMD ["-h"]
