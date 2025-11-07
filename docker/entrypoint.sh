#!/bin/sh
set -e

if [ -n "$POSTGRES_HOST" ]; then
  until pg_isready -h "$POSTGRES_HOST" -p "${POSTGRES_PORT:-5432}" -U "${POSTGRES_USER:-postgres}" >/dev/null 2>&1; do
    echo "Waiting for PostgreSQL to become available..."
    sleep 1
  done
fi

# Check if database has existing tables (from backup restore)
if python manage.py showmigrations --plan | grep -q "^\[X\]"; then
    echo "Database appears to have existing migrations, skipping migration..."
else
    # Try normal migration first
    if ! python manage.py migrate --noinput; then
        echo "Migration failed, trying --fake-initial for restored database..."
        python manage.py migrate --fake-initial --noinput
    fi
fi

python manage.py collectstatic --noinput

exec "$@"