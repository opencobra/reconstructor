"""
MATLAB HTTP Client

This module provides a client to communicate with the MATLAB HTTP API server
running in a separate Docker container.
"""

import os
import requests
import logging
from typing import Dict, Any, List, Optional

logger = logging.getLogger(__name__)


class MatlabHTTPClient:
    """
    Client to communicate with MATLAB HTTP API server.
    Replaces direct MATLAB Engine API usage with HTTP requests.
    """
    _instance: Optional['MatlabHTTPClient'] = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(MatlabHTTPClient, cls).__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        """Initialize the HTTP client."""
        if self._initialized:
            return
        
        # Get MATLAB server configuration from environment
        self.matlab_host = os.getenv('MATLAB_HOST', 'matlab')
        self.matlab_port = int(os.getenv('MATLAB_HTTP_PORT', '9090'))
        self.base_url = f"http://{self.matlab_host}:{self.matlab_port}"
        self.timeout = 300  # 5 minutes timeout for MATLAB operations
        
        logger.info(f"MATLAB HTTP Client initialized: {self.base_url}")
        self._initialized = True

    def health_check(self) -> bool:
        """
        Check if the MATLAB server is healthy and responsive.
        
        Returns:
            bool: True if server is healthy, False otherwise
        """
        try:
            response = requests.get(
                f"{self.base_url}/health",
                timeout=5
            )
            return response.status_code == 200
        except Exception as e:
            logger.error(f"MATLAB health check failed: {e}")
            return False

    def execute(
        self,
        function: str,
        *args,
        nargout: int = 1,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Execute a MATLAB function.
        
        Args:
            function: Name of the MATLAB function to execute
            *args: Positional arguments to pass to the function
            nargout: Number of output arguments (default 1)
            **kwargs: Keyword arguments to pass to the function
            
        Returns:
            Dictionary with 'status' and either 'result' or 'message':
            - {'status': 'success', 'result': <return_value>}
            - {'status': 'error', 'message': <error_message>}
        """
        try:
            payload = {
                'function': function,
                'args': list(args),
                'kwargs': kwargs,
                'nargout': nargout
            }
            
            logger.debug(f"Executing MATLAB function: {function}")
            
            response = requests.post(
                f"{self.base_url}/execute",
                json=payload,
                timeout=self.timeout
            )
            
            result = response.json()
            
            if response.status_code == 200:
                logger.debug(f"Function executed successfully: {function}")
                return result
            else:
                logger.error(f"Function execution failed: {result.get('message', 'Unknown error')}")
                return result
                
        except requests.exceptions.Timeout:
            error_msg = f"MATLAB function {function} timed out after {self.timeout} seconds"
            logger.error(error_msg)
            return {'status': 'error', 'message': error_msg}
        except requests.exceptions.ConnectionError as e:
            error_msg = f"Cannot connect to MATLAB server at {self.base_url}: {str(e)}"
            logger.error(error_msg)
            return {'status': 'error', 'message': error_msg}
        except Exception as e:
            error_msg = f"Unexpected error executing {function}: {str(e)}"
            logger.error(error_msg)
            return {'status': 'error', 'message': error_msg}

    def eval(self, code: str, nargout: int = 0) -> Dict[str, Any]:
        """
        Evaluate MATLAB code.
        
        Args:
            code: MATLAB code string to evaluate
            nargout: Number of output arguments (default 0)
            
        Returns:
            Dictionary with 'status' and either 'result' or 'message':
            - {'status': 'success', 'result': <return_value>}  (if nargout > 0)
            - {'status': 'success'}  (if nargout == 0)
            - {'status': 'error', 'message': <error_message>}
        """
        try:
            payload = {
                'code': code,
                'nargout': nargout
            }
            
            logger.debug(f"Evaluating MATLAB code: {code[:100]}...")
            
            response = requests.post(
                f"{self.base_url}/eval",
                json=payload,
                timeout=self.timeout
            )
            
            result = response.json()
            
            if response.status_code == 200:
                logger.debug("Code evaluated successfully")
                return result
            else:
                logger.error(f"Code evaluation failed: {result.get('message', 'Unknown error')}")
                return result
                
        except requests.exceptions.Timeout:
            error_msg = f"MATLAB code evaluation timed out after {self.timeout} seconds"
            logger.error(error_msg)
            return {'status': 'error', 'message': error_msg}
        except requests.exceptions.ConnectionError as e:
            error_msg = f"Cannot connect to MATLAB server at {self.base_url}: {str(e)}"
            logger.error(error_msg)
            return {'status': 'error', 'message': error_msg}
        except Exception as e:
            error_msg = f"Unexpected error evaluating code: {str(e)}"
            logger.error(error_msg)
            return {'status': 'error', 'message': error_msg}

    @property
    def is_connected(self) -> bool:
        """Check if the MATLAB server is connected and responsive."""
        return self.health_check()


# Alias for backwards compatibility
MatlabSessionManager = MatlabHTTPClient
