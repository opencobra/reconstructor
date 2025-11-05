#!/bin/sh
set -e

if [ -n "$POSTGRES_HOST" ]; then
  until pg_isready -h "$POSTGRES_HOST" -p "${POSTGRES_PORT:-5432}" -U "${POSTGRES_USER:-postgres}" >/dev/null 2>&1; do
    echo "Waiting for PostgreSQL to become available..."
    sleep 1
  done
fi

python manage.py migrate --noinput
python manage.py collectstatic --noinput

exec "$@"
