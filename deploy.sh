#/home/saleh/reconstructor/.venv/bin/python
set -euo pipefail
cd reconstructor  

echo "Pulling latest code…"
git fetch --prune
git checkout main
git reset --hard origin/main

echo "Installing deps…"
source .venv//bin/activate
pip install -r requirements.txt

echo "Running migrations & collectstatic…"
cd curationTool/
python manage.py migrate --noinput
python manage.py collectstatic --noinput

echo "Reloading services…"
sudo systemctl daemon-reload
sudo systemctl restart gunicorn.socket gunicorn.service
sudo nginx -t && sudo systemctl reload nginx

echo "✓ Deploy finished at $(date -Iseconds)"