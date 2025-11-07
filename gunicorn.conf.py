"""
Gunicorn configuration file for production Django deployment.
This configuration ensures all errors and logs are visible in Docker logs.
"""

import multiprocessing

# Server socket
bind = "0.0.0.0:8000"
backlog = 2048

# Worker processes
workers = multiprocessing.cpu_count() * 2 + 1
worker_class = "sync"
worker_connections = 1000
timeout = 120
keepalive = 2

# Logging
# Access log - set to '-' to log to stdout
accesslog = "-"
# Error log - set to '-' to log to stderr
errorlog = "-"
# Log level - set to 'debug' for maximum verbosity, 'info' for production
loglevel = "info"
# Capture stdout/stderr from application
capture_output = True
# Enable detailed error logging
enable_stdio_inheritance = True

# Process naming
proc_name = "curationTool"

# Server mechanics
daemon = False
pidfile = None
umask = 0
user = None
group = None
tmp_upload_dir = None

# Logging format
access_log_format = '%(h)s %(l)s %(u)s %(t)s "%(r)s" %(s)s %(b)s "%(f)s" "%(a)s"'
