# syntax=docker/dockerfile:1

# Stage 1: Install MATLAB Engine API in MATLAB container
# The mathworks/matlab:r2024b image uses Python 3.12
FROM mathworks/matlab:r2024b AS matlab-builder
USER root
RUN apt-get update && apt-get install -y python3 python3-pip python3-dev && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/matlab/R2024b/extern/engines/python
# Install the MATLAB Engine API into a staging directory that we can copy later
RUN python3 setup.py install --prefix=/tmp/matlabengine
RUN find /tmp/matlabengine -maxdepth 4 -type d

# Stage 2: Web application
# Must match the Python version used in matlab-builder (3.12)
FROM python:3.12-slim AS base

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

# Copy MATLAB Engine API from builder - installed under local/ prefix
# The path structure is /tmp/matlabengine/local/lib/python3.12/dist-packages/
COPY --from=matlab-builder /tmp/matlabengine/local/ /usr/local/

# Copy application code
COPY . .

RUN chmod +x docker/entrypoint.sh

WORKDIR /app/curationTool
ENV DJANGO_SETTINGS_MODULE=reactions_project.settings

EXPOSE 8000

ENTRYPOINT ["/app/docker/entrypoint.sh"]
CMD ["gunicorn", "reactions_project.wsgi:application", "--config", "/app/gunicorn.conf.py"]
