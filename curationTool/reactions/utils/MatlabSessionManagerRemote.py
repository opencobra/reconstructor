"""
Remote MATLAB Session Manager

This module provides a singleton class to manage connections to a remote MATLAB
Engine API server running in a separate Docker container. It replaces the local
MATLAB session with a network-based connection.

Usage:
    from reactions.utils.MatlabSessionManagerRemote import MatlabSessionManager
    
    matlab_session = MatlabSessionManager()
    result = matlab_session.execute('functionName', arg1, arg2, kwarg1=value1)
"""

import os
import time
import logging
import matlab.engine
from typing import Dict, Any, Optional

logger = logging.getLogger(__name__)


class MatlabSessionManager:
    """
    Singleton class to manage a remote MATLAB Engine session.
    
    This connects to a MATLAB Engine API server running in a separate Docker container
    instead of starting a local MATLAB instance.
    """
    _instance: Optional['MatlabSessionManager'] = None
    _max_retries: int = 3
    _retry_delay: int = 5  # seconds

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(MatlabSessionManager, cls).__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        """Initialize the remote MATLAB session connection."""
        if self._initialized:
            return
            
        self.engine = None
        self.session_name = os.getenv('MATLAB_SESSION_NAME', 'matlab_shared_session')
        self.use_remote = os.getenv('MATLAB_REMOTE_ENABLED', 'false').lower() == 'true'
        
        if self.use_remote:
            self._connect_to_remote_matlab()
        else:
            # Fallback to local MATLAB (original behavior)
            self._connect_to_local_matlab()
        
        self._initialized = True

    def _connect_to_local_matlab(self):
        """
        Connect to local MATLAB instance (legacy mode).
        This is the original implementation for backwards compatibility.
        """
        try:
            logger.info("Starting local MATLAB session...")
            self.engine = matlab.engine.start_matlab()
            self._setup_cobra_toolbox()
            logger.info("Local MATLAB session started successfully")
        except Exception as e:
            logger.error(f"Failed to start local MATLAB session: {e}")
            raise

    def _connect_to_remote_matlab(self):
        """
        Connect to remote MATLAB Engine API server running in Docker container.
        Implements retry logic for connection reliability.
        """
        for attempt in range(1, self._max_retries + 1):
            try:
                logger.info(f"Attempting to connect to remote MATLAB session '{self.session_name}' (attempt {attempt}/{self._max_retries})...")
                
                # Find shared MATLAB sessions
                available_sessions = matlab.engine.find_matlab()
                logger.info(f"Available MATLAB sessions: {available_sessions}")
                
                if self.session_name not in available_sessions:
                    raise ConnectionError(
                        f"MATLAB session '{self.session_name}' not found. "
                        f"Available sessions: {available_sessions}"
                    )
                
                # Connect to the shared session
                self.engine = matlab.engine.connect_matlab(self.session_name)
                logger.info(f"Successfully connected to remote MATLAB session: {self.session_name}")
                
                # Test the connection
                test_result = self.engine.eval("1+1", nargout=1)
                logger.info(f"Connection test successful: 1+1={test_result}")
                
                return  # Success
                
            except Exception as e:
                logger.warning(f"Connection attempt {attempt} failed: {e}")
                
                if attempt < self._max_retries:
                    logger.info(f"Retrying in {self._retry_delay} seconds...")
                    time.sleep(self._retry_delay)
                else:
                    logger.error(f"Failed to connect to remote MATLAB after {self._max_retries} attempts")
                    raise ConnectionError(
                        f"Unable to connect to remote MATLAB session '{self.session_name}'. "
                        "Ensure the MATLAB container is running and the session is shared."
                    ) from e

    def _setup_cobra_toolbox(self):
        """
        Set up COBRA Toolbox and custom script directories.
        Only used in local mode; remote mode handles this on server startup.
        """
        if self.use_remote:
            logger.info("Skipping COBRA setup (handled by remote MATLAB server)")
            return
            
        # Read configuration from environment variables (only for local mode)
        script_directories_env = os.getenv('SCRIPT_DIRECTORIES', '')
        script_directories = [d.strip() for d in script_directories_env.split(',') if d.strip()]
        cobra_path = os.getenv('COBRA_PATH', '')

        for script_directory in script_directories:
            logger.info(f"Adding MATLAB path: {script_directory}")
            self.engine.addpath(script_directory, nargout=0)
            
        if cobra_path:
            logger.info(f"Adding COBRA Toolbox path: {cobra_path}")
            self.engine.addpath(cobra_path, nargout=0)
            logger.info("Initializing COBRA Toolbox...")
            self.engine.eval("initCobraToolbox(0)", nargout=0)

    def execute(self, command: str, *args, **kwargs) -> Dict[str, Any]:
        """
        Execute a MATLAB command or function.
        
        Args:
            command: Name of the MATLAB function to execute
            *args: Positional arguments to pass to the function
            **kwargs: Keyword arguments to pass to the function
            
        Returns:
            Dictionary with 'status' and either 'result' or 'message':
            - {'status': 'success', 'result': <return_value>}
            - {'status': 'error', 'message': <error_message>}
        """
        if not self.engine:
            return {
                'status': 'error',
                'message': 'MATLAB engine not initialized'
            }
            
        try:
            if hasattr(self.engine, command):
                logger.debug(f"Executing MATLAB command: {command}")
                matlab_function = getattr(self.engine, command)
                result = matlab_function(*args, **kwargs)
                logger.debug(f"Command executed successfully: {command}")
                return {'status': 'success', 'result': result}
            else:
                error_msg = f'Command {command} not found in MATLAB engine'
                logger.error(error_msg)
                return {'status': 'error', 'message': error_msg}
                
        except matlab.engine.MatlabExecutionError as e:
            error_msg = f"MATLAB execution error in {command}: {str(e)}"
            logger.error(error_msg)
            return {'status': 'error', 'message': str(e)}
        except Exception as e:
            error_msg = f"Unexpected error executing {command}: {str(e)}"
            logger.error(error_msg)
            return {'status': 'error', 'message': str(e)}

    def quit(self):
        """
        Disconnect from the MATLAB session.
        Note: For remote sessions, this only disconnects the client;
        the server continues running.
        """
        if self.engine:
            try:
                if self.use_remote:
                    logger.info("Disconnecting from remote MATLAB session...")
                    # For remote connections, we just disconnect, don't quit the server
                    self.engine = None
                else:
                    logger.info("Quitting local MATLAB session...")
                    self.engine.quit()
                    self.__class__._instance = None
            except Exception as e:
                logger.error(f"Error while disconnecting from MATLAB: {e}")

    def reconnect(self):
        """
        Reconnect to the MATLAB session.
        Useful if the connection is lost.
        """
        logger.info("Attempting to reconnect to MATLAB...")
        self.quit()
        self.__class__._instance = None
        self._initialized = False
        self.__init__()

    @property
    def is_connected(self) -> bool:
        """Check if the MATLAB engine is connected and responsive."""
        if not self.engine:
            return False
        try:
            self.engine.eval("1+1", nargout=0)
            return True
        except Exception:
            return False
