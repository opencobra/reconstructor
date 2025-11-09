# MATLAB Container with MATLAB Engine API
FROM mathworks/matlab:r2024b

# Install Python and necessary packages for MATLAB Engine API
USER root

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install MATLAB Engine API for Python
# The MATLAB Engine API is located in the MATLAB installation
WORKDIR /opt/matlab/R2024b/extern/engines/python
RUN python3 setup.py install

# Create directories for MATLAB scripts and toolboxes
RUN mkdir -p /matlab/scripts
RUN mkdir -p /matlab/toolboxes/cobratoolbox
RUN mkdir -p /matlab/toolboxes/vmh_revamped

# Create MATLAB user directories to prevent startup errors
RUN mkdir -p /root/Documents/MATLAB

WORKDIR /matlab/scripts

# Copy any custom MATLAB scripts or toolboxes
# COPY ./scripts /matlab/scripts

# Create a startup script that runs MATLAB in server mode
COPY matlab_engine_server.py /matlab/
RUN chmod +x /matlab/matlab_engine_server.py

# Create directory for license files
RUN mkdir -p /licenses

# Expose port for MATLAB Engine API
EXPOSE 9090

# Stay as root user since license is activated for root
USER root

# Override the default MATLAB entrypoint to run our Python script instead
ENTRYPOINT []

# Run MATLAB Engine API server
CMD ["python3", "/matlab/matlab_engine_server.py"]