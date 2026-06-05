# Docker deployment

## Postgres data

The Postgres container stores its database files in a host bind mount:

```yaml
${POSTGRES_DATA_DIR:-./docker/postgres_data}:/var/lib/postgresql/data
```

That means rebuilding images and recreating containers does not reset the Django
database. The data lives on the server filesystem and the Postgres container uses
that directory.

For production, set an absolute host path in `.env`:

```env
POSTGRES_DATA_DIR=/home/saleh/reconstructor/docker/postgres_data
```

You can also move it outside the repo, for example:

```env
POSTGRES_DATA_DIR=/srv/reconstructor/postgres_data
```

If you move the directory, stop Compose first, copy the existing
`docker/postgres_data` directory to the new path while preserving ownership and
permissions, update `.env`, then start Compose again.

Do not remove the data directory. `docker compose build`, `docker compose up -d`,
and `docker compose down` should not delete it. Avoid destructive cleanup commands
such as `rm -rf docker/postgres_data`.

## Deploy sequence

After pulling new code:

```bash
sudo docker compose build
sudo docker compose up -d
```

The web container entrypoint runs these before Gunicorn starts:

```bash
python manage.py migrate --noinput
python manage.py collectstatic --noinput
```

So pending Django migrations and static files are applied automatically during
`up -d` whenever the web container starts or is recreated.

## Useful checks

```bash
sudo docker compose ps
sudo docker compose logs -f web
sudo docker compose logs -f db
```

To confirm migrations are applied:

```bash
sudo docker compose exec web python manage.py showmigrations reactions
```
