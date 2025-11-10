# syntax=docker/dockerfile:1

# Stage 1: Install MATLAB Engine API in MATLAB container
# The mathworks/matlab:r2024b image uses Python 3.12 (Ubuntu Noble)
FROM mathworks/matlab:r2024b AS matlab-builder
USER root
RUN apt-get update && apt-get install -y python3 python3-pip python3-dev && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/matlab/R2024b/extern/engines/python
# Build with system Python 3.12; the compiled extension uses stable ABI and works with 3.11+
RUN python3 -m pip install --upgrade pip setuptools wheel && \
    python3 setup.py install --prefix=/tmp/matlabengine
RUN find /tmp/matlabengine -maxdepth 5 -type d -print
# The installer creates lib/python3.12/site-packages but we'll copy to 3.11 runtime
# Create a symlink so both 3.11 and 3.12 paths work
RUN if [ -d /tmp/matlabengine/lib/python3.12 ]; then \
        cd /tmp/matlabengine/lib && \
        ln -s python3.12 python3.11; \
    fi

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
