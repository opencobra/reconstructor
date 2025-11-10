# syntax=docker/dockerfile:1

# Stage 1: Install MATLAB Engine API in MATLAB container
FROM mathworks/matlab:r2024b AS matlab-builder
USER root
RUN apt-get update && apt-get install -y python3 python3-pip python3-dev && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/matlab/R2024b/extern/engines/python
# Install the MATLAB Engine API (this installs all dependencies including compiled extensions)
RUN python3 setup.py install

# Stage 2: Web application
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
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy MATLAB Engine API and all its dependencies from builder stage
# This includes matlab, matlabengineforpython, and matlabmultidimarrayforpython packages
COPY --from=matlab-builder /usr/local/lib/python3.11/site-packages/matlab /usr/local/lib/python3.11/site-packages/matlab
COPY --from=matlab-builder /usr/local/lib/python3.11/site-packages/matlabengineforpython* /usr/local/lib/python3.11/site-packages/
COPY --from=matlab-builder /usr/local/lib/python3.11/site-packages/matlabmultidimarrayforpython* /usr/local/lib/python3.11/site-packages/

# Copy application code
COPY . .

RUN chmod +x docker/entrypoint.sh

WORKDIR /app/curationTool
ENV DJANGO_SETTINGS_MODULE=reactions_project.settings

EXPOSE 8000

ENTRYPOINT ["/app/docker/entrypoint.sh"]
CMD ["gunicorn", "reactions_project.wsgi:application", "--config", "/app/gunicorn.conf.py"]
