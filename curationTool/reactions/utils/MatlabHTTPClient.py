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
        
        print(f"[INFO MatlabHTTPClient] Initializing with host={self.matlab_host}, port={self.matlab_port}", flush=True)
        print(f"[INFO MatlabHTTPClient] Base URL: {self.base_url}", flush=True)
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

    def execute(self, function_name: str, *args, **kwargs) -> Dict[str, Any]:
        """
        Execute a MATLAB function with arguments.
        
        Args:
            function_name: Name of the MATLAB function to execute
            *args: Positional arguments for the function
            **kwargs: Keyword arguments for the function
            
        Returns:
            Dict with 'status' ('success' or 'error'), 'result', and optional 'message'
        """
        print(f"[DEBUG MatlabHTTPClient] Executing MATLAB function: {function_name}", flush=True)
        print(f"[DEBUG MatlabHTTPClient] Args: {args}, Kwargs: {kwargs}", flush=True)
        print(f"[DEBUG MatlabHTTPClient] Sending request to: {self.base_url}/execute", flush=True)
        
        try:
            response = requests.post(
                f"{self.base_url}/execute",
                json={
                    'function': function_name,
                    'args': list(args),
                    'kwargs': kwargs
                },
                timeout=self.timeout
            )
            print(f"[DEBUG MatlabHTTPClient] Response status: {response.status_code}", flush=True)
            print(f"[DEBUG MatlabHTTPClient] Response body: {response.text[:200]}", flush=True)
            
            response.raise_for_status()
            result = response.json()
            print(f"[DEBUG MatlabHTTPClient] Result: {result}", flush=True)
            return result
            
        except requests.exceptions.ConnectionError as e:
            error_msg = f"Cannot connect to MATLAB server at {self.base_url}: {e}"
            print(f"[ERROR MatlabHTTPClient] {error_msg}", flush=True)
            logger.error(error_msg)
            return {'status': 'error', 'message': error_msg}
        except requests.exceptions.Timeout as e:
            error_msg = f"Request to MATLAB server timed out after {self.timeout}s: {e}"
            print(f"[ERROR MatlabHTTPClient] {error_msg}", flush=True)
            logger.error(error_msg)
            return {'status': 'error', 'message': error_msg}
        except Exception as e:
            error_msg = f"Error executing MATLAB function '{function_name}': {e}"
            print(f"[ERROR MatlabHTTPClient] {error_msg}", flush=True)
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
