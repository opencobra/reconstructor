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
timeout = 4500
keepalive = 2

# Logging
# Access log - set to '-' to log to stdout
accesslog = "-"
# Error log - set to '-' to log to stderr
errorlog = "-"
# Log level - set to 'debug' for maximum verbosity, 'info' for production
loglevel = "debug"
# Capture stdout/stderr from application
capture_output = True
# Enable detailed error logging
enable_stdio_inheritance = True
# Log to stdout/stderr instead of files
logconfig_dict = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'generic': {
            'format': '[%(levelname)s] %(asctime)s [%(process)d] [%(name)s] %(message)s',
            'datefmt': '%Y-%m-%d %H:%M:%S',
        },
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'generic',
            'stream': 'ext://sys.stdout',
        },
        'error_console': {
            'class': 'logging.StreamHandler',
            'formatter': 'generic',
            'stream': 'ext://sys.stderr',
        },
    },
    'root': {
        'level': 'INFO',
        'handlers': ['console'],
    },
    'loggers': {
        'gunicorn.error': {
            'level': 'DEBUG',
            'handlers': ['error_console'],
            'propagate': False,
        },
        'gunicorn.access': {
            'level': 'INFO',
            'handlers': ['console'],
            'propagate': False,
        },
    },
}


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
