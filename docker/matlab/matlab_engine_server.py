#!/usr/bin/env python3
"""
MATLAB Engine API Server
This script starts a MATLAB session and keeps it running,
making it available for remote connections via the MATLAB Engine API.
"""
import os
import sys
import time
import matlab.engine
from pathlib import Path

def start_matlab_engine_server():
    """
    Start MATLAB Engine in shared mode so it can be connected to remotely.
    """
    print("Starting MATLAB Engine Server...", flush=True)
    
    try:
        # Start MATLAB with a specific session name for remote connections
        engine = matlab.engine.start_matlab()
        
        # Setup COBRA Toolbox if path is provided
        cobra_path = os.getenv('COBRA_PATH', '')
        script_directories_env = os.getenv('SCRIPT_DIRECTORIES', '')
        script_directories = [d.strip() for d in script_directories_env.split(',') if d.strip()]
        
        print("Setting up MATLAB environment...", flush=True)
        
        # Add script directories to MATLAB path
        for script_directory in script_directories:
            if script_directory and Path(script_directory).exists():
                print(f"Adding path: {script_directory}", flush=True)
                engine.addpath(script_directory, nargout=0)
        
        # Initialize COBRA Toolbox if available
        if cobra_path and Path(cobra_path).exists():
            print(f"Adding COBRA path: {cobra_path}", flush=True)
            engine.addpath(cobra_path, nargout=0)
            print("Initializing COBRA Toolbox...", flush=True)
            engine.eval("initCobraToolbox(0)", nargout=0)
        
        # Share the MATLAB session with a specific name
        session_name = os.getenv('MATLAB_SESSION_NAME', 'matlab_shared_session')
        print(f"Sharing MATLAB session as: {session_name}", flush=True)
        engine.matlab.engine.shareEngine(session_name, nargout=0)
        
        print("MATLAB Engine Server is ready and accepting connections!", flush=True)
        print(f"Session name: {session_name}", flush=True)
        
        # Keep the server running
        while True:
            time.sleep(10)
            # Optionally perform health checks
            try:
                engine.eval("1+1", nargout=0)
            except Exception as e:
                print(f"Health check failed: {e}", flush=True)
                break
                
    except Exception as e:
        print(f"Failed to start MATLAB Engine Server: {e}", file=sys.stderr, flush=True)
        sys.exit(1)

if __name__ == "__main__":
    start_matlab_engine_server()
