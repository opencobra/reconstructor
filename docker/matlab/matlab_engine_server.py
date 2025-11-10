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
            try:
                # Initialize without checking for updates (since we're in a container)
                engine.eval("initCobraToolbox(false)", nargout=0)
                print("COBRA Toolbox initialized successfully", flush=True)
            except Exception as e:
                print(f"Warning: COBRA Toolbox initialization failed: {e}", flush=True)
                print("Continuing without COBRA Toolbox...", flush=True)
        
        # Share the MATLAB session with a specific name
        session_name = os.getenv('MATLAB_SESSION_NAME', 'matlab_shared_session')
        print(f"Sharing MATLAB session as: {session_name}", flush=True)
        
        try:
            engine.matlab.engine.shareEngine(session_name, nargout=0)
            print("MATLAB Engine Server is ready and accepting connections!", flush=True)
            print(f"Session name: {session_name}", flush=True)
        except Exception as share_error:
            # If session name already exists, it might be from a previous container run
            print(f"Warning: Could not share session as '{session_name}': {share_error}", flush=True)
            
            # Check if we can find and connect to the existing session
            try:
                existing_sessions = matlab.engine.find_matlab()
                print(f"Found existing MATLAB sessions: {existing_sessions}", flush=True)
                
                if session_name in existing_sessions:
                    print(f"Session '{session_name}' already exists. Attempting to verify it's alive...", flush=True)
                    try:
                        # Try to connect to the existing session to see if it's valid
                        test_engine = matlab.engine.connect_matlab(session_name)
                        test_engine.eval("1+1", nargout=0)
                        test_engine.quit()
                        print(f"Existing session '{session_name}' is alive and functional.", flush=True)
                        print("The existing session will continue to serve requests.", flush=True)
                        print("This container will exit to avoid conflicts.", flush=True)
                        sys.exit(0)  # Exit successfully since a valid session exists
                    except Exception as connect_error:
                        print(f"Existing session appears to be dead/orphaned: {connect_error}", flush=True)
                        print("Waiting 5 seconds for session cleanup and retrying...", flush=True)
                        time.sleep(5)
                        
                        # Retry sharing after waiting
                        try:
                            engine.matlab.engine.shareEngine(session_name, nargout=0)
                            print(f"Successfully shared session as '{session_name}' after retry", flush=True)
                        except Exception as retry_error:
                            print(f"Retry failed: {retry_error}. Using unique session name instead.", flush=True)
                            import uuid
                            fallback_name = f"{session_name}_{uuid.uuid4().hex[:8]}"
                            engine.matlab.engine.shareEngine(fallback_name, nargout=0)
                            print(f"Using fallback session name: {fallback_name}", flush=True)
                else:
                    # Session name not in list, but shareEngine still failed - use fallback
                    import uuid
                    fallback_name = f"{session_name}_{uuid.uuid4().hex[:8]}"
                    print(f"Sharing with fallback name: {fallback_name}", flush=True)
                    engine.matlab.engine.shareEngine(fallback_name, nargout=0)
                    print(f"MATLAB Engine Server is ready with session name: {fallback_name}", flush=True)
            except Exception as find_error:
                print(f"Error checking existing sessions: {find_error}", flush=True)
                print("Continuing with unnamed shared session...", flush=True)
        
        # Keep the server running
        while True:
            time.sleep(10)
            # Perform health check (suppress output with nargout=0 and semicolon)
            try:
                # Use nargout=0 and eval with semicolon to suppress output
                engine.eval("x = 1 + 1;", nargout=0)
            except Exception as e:
                print(f"Health check failed: {e}", flush=True)
                break
                
    except Exception as e:
        print(f"Failed to start MATLAB Engine Server: {e}", file=sys.stderr, flush=True)
        sys.exit(1)

if __name__ == "__main__":
    start_matlab_engine_server()
