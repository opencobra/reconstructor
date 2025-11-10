# syntax=docker/dockerfile:1

# Stage 1: Install MATLAB Engine API in MATLAB container
# The mathworks/matlab:r2024b image uses Python 3.12
FROM mathworks/matlab:r2024b AS matlab-builder
USER root
# Install Python 3.11 in the MATLAB image and use it to build the MATLAB Engine
# so the engine wheel/extension matches the Python runtime we use for the web image.
RUN apt-get update && apt-get install -y \
    python3.11 \
    python3.11-dev \
    python3.11-venv \
    python3-pip \
    python3-distutils \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/matlab/R2024b/extern/engines/python
# Ensure pip/setuptools/wheel for python3.11 are available, then build the engine
RUN python3.11 -m pip install --upgrade pip setuptools wheel && \
    python3.11 setup.py install --prefix=/tmp/matlabengine
# Inspect created layout
RUN find /tmp/matlabengine -maxdepth 5 -type d -print

# Stage 2: Web application
# We'll use Python 3.11 for the web image so manylinux wheels for your requirements
FROM python:3.11-slim AS base

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libpq-dev \
    default-jre \
    postgresql-client \
    libgl1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python dependencies first (for better caching)
COPY requirements.txt ./
# Upgrade pip/setuptools/wheel first so build backends can be imported if needed
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt

# Copy MATLAB Engine API from builder - installed under local/ prefix
# The path structure is /tmp/matlabengine/local/lib/python3.12/dist-packages/
COPY --from=matlab-builder /tmp/matlabengine/ /usr/local/

# Copy application code
COPY . .

RUN chmod +x docker/entrypoint.sh

WORKDIR /app/curationTool
ENV DJANGO_SETTINGS_MODULE=reactions_project.settings

EXPOSE 8000

ENTRYPOINT ["/app/docker/entrypoint.sh"]
CMD ["gunicorn", "reactions_project.wsgi:application", "--config", "/app/gunicorn.conf.py"]
