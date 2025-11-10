# syntax=docker/dockerfile:1

# Stage 1: Install MATLAB Engine API in MATLAB container
# Keep PYTHON_VERSION aligned with the interpreter used to build the MATLAB Engine
# so the compiled extension matches the runtime ABI in the web image.
ARG PYTHON_VERSION=3.10
FROM mathworks/matlab:r2024b AS matlab-builder
USER root
RUN apt-get update && apt-get install -y python3 python3-pip python3-dev && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/matlab/R2024b/extern/engines/python
# Install the MATLAB Engine API into a staging directory that we can copy later
RUN python3 setup.py install --prefix=/tmp/matlabengine
RUN find /tmp/matlabengine -maxdepth 4 -type d

# Stage 2: Web application
FROM python:${PYTHON_VERSION}-slim AS base
ARG PYTHON_VERSION

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
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy MATLAB Engine API and all its dependencies from builder stage
# These are installed under the staging prefix
COPY --from=matlab-builder /tmp/matlabengine/lib/python${PYTHON_VERSION}/site-packages/ /usr/local/lib/python${PYTHON_VERSION}/site-packages/

# Copy application code
COPY . .

RUN chmod +x docker/entrypoint.sh

WORKDIR /app/curationTool
ENV DJANGO_SETTINGS_MODULE=reactions_project.settings

EXPOSE 8000

ENTRYPOINT ["/app/docker/entrypoint.sh"]
CMD ["gunicorn", "reactions_project.wsgi:application", "--config", "/app/gunicorn.conf.py"]
