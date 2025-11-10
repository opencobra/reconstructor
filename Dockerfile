# syntax=docker/dockerfile:1

# Web application - no MATLAB Engine needed
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

# Copy application code
COPY . .

RUN chmod +x docker/entrypoint.sh

WORKDIR /app/curationTool
ENV DJANGO_SETTINGS_MODULE=reactions_project.settings

EXPOSE 8000

ENTRYPOINT ["/app/docker/entrypoint.sh"]
CMD ["gunicorn", "reactions_project.wsgi:application", "--config", "/app/gunicorn.conf.py"]
