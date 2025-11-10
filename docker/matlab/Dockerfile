# MATLAB Container with MATLAB Engine API
FROM mathworks/matlab:r2024b

# Install Python and necessary packages for MATLAB Engine API
USER root

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install MATLAB Engine API for Python
# The MATLAB Engine API is located in the MATLAB installation
WORKDIR /opt/matlab/R2024b/extern/engines/python
RUN python3 setup.py install

# Install Flask for HTTP API server
RUN python3 -m pip install --break-system-packages flask

# Create directories for MATLAB scripts and toolboxes
RUN mkdir -p /matlab/scripts
RUN mkdir -p /matlab/toolboxes/vmh_revamped
RUN mkdir -p /matlab/data

# Clone and initialize COBRA Toolbox with submodules
RUN git clone --depth=1 https://github.com/opencobra/cobratoolbox.git /matlab/toolboxes/cobratoolbox && \
    cd /matlab/toolboxes/cobratoolbox && \
    git submodule update --init --recursive

# Create MATLAB user directories to prevent startup errors
RUN mkdir -p /root/Documents/MATLAB

# Make /matlab directory writable for MATLAB operations
RUN chmod -R 777 /matlab

WORKDIR /matlab/scripts

# Note: vmh_revamped should be accessible in build context
# If using symlink: ln -s /home/saleh/vmh_revamped vmh_revamped (from project root)
# Copy VMH revamped code from build context
COPY --chown=root:root vmh_revamped/ /matlab/toolboxes/vmh_revamped/

# Copy any custom MATLAB scripts or toolboxes  
# COPY ./scripts /matlab/scripts

# Create a startup script that runs MATLAB in server mode
COPY docker/matlab/matlab_http_server.py /matlab/
RUN chmod +x /matlab/matlab_http_server.py

# Create directory for license files
RUN mkdir -p /licenses

# Expose port for MATLAB Engine API
EXPOSE 9090

# Stay as root user since license is activated for root
USER root

# Override the default MATLAB entrypoint to run our Python script instead
ENTRYPOINT []

# Run MATLAB HTTP API server
CMD ["python3", "/matlab/matlab_http_server.py"]