#!/usr/bin/env python3
"""
MATLAB HTTP API Server
Provides a REST API to execute MATLAB functions remotely without requiring
the MATLAB Engine API in client containers.
"""
import os
import sys
import atexit
import traceback
from pathlib import Path
from flask import Flask, request, jsonify
import matlab.engine

app = Flask(__name__)
matlab_engine = None

def initialize_matlab():
    """Initialize MATLAB engine and set up paths"""
    global matlab_engine
    
    print("Starting MATLAB Engine...", flush=True)
    matlab_engine = matlab.engine.start_matlab()
    
    # Setup environment
    cobra_path = os.getenv('COBRA_PATH', '')
    script_directories_env = os.getenv('SCRIPT_DIRECTORIES', '')
    script_directories = [d.strip() for d in script_directories_env.split(',') if d.strip()]
    
    print("Setting up MATLAB environment...", flush=True)
    
    # Add script directories
    for script_directory in script_directories:
        if script_directory and Path(script_directory).exists():
            print(f"Adding path: {script_directory}", flush=True)
            matlab_engine.addpath(script_directory, nargout=0)
    
    # Initialize COBRA Toolbox
    if cobra_path and Path(cobra_path).exists():
        print(f"Adding COBRA path: {cobra_path}", flush=True)
        matlab_engine.addpath(cobra_path, nargout=0)
        
        print("Initializing COBRA Toolbox...", flush=True)
        try:
            matlab_engine.eval("initCobraToolbox(false)", nargout=0)
            print("COBRA Toolbox initialized successfully", flush=True)
        except Exception as e:
            print(f"Warning: COBRA Toolbox initialization failed: {e}", flush=True)
    
    print("MATLAB Engine initialized and ready!", flush=True)

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    try:
        if matlab_engine is None:
            return jsonify({'status': 'error', 'message': 'MATLAB engine not initialized'}), 503
        
        # Quick health check
        matlab_engine.eval("1+1", nargout=0)
        return jsonify({'status': 'ok', 'message': 'MATLAB engine is healthy'}), 200
    except Exception as e:
        return jsonify({'status': 'error', 'message': str(e)}), 500

@app.route('/execute', methods=['POST'])
def execute_function():
    """
    Execute a MATLAB function
    
    Request body:
    {
        "function": "functionName",
        "args": [arg1, arg2, ...],  // optional positional arguments
        "kwargs": {key1: val1, ...},  // optional keyword arguments
        "nargout": 1  // optional, number of output arguments (default 1)
    }
    
    Response:
    {
        "status": "success" | "error",
        "result": <result_value>,  // if status is success
        "message": <error_message>  // if status is error
    }
    """
    try:
        data = request.get_json()
        
        if not data or 'function' not in data:
            return jsonify({
                'status': 'error',
                'message': 'Missing required field: function'
            }), 400
        
        function_name = data['function']
        args = data.get('args', [])
        kwargs = data.get('kwargs', {})
        nargout = data.get('nargout', 1)

        print(
            f"[MATLAB HTTP] Received execute request: function={function_name}, args={args}, kwargs={kwargs}, nargout={nargout}",
            flush=True
        )
        
        # Validate that function exists
        if not hasattr(matlab_engine, function_name):
            return jsonify({
                'status': 'error',
                'message': f'Function {function_name} not found in MATLAB engine'
            }), 404
        
        # Get the MATLAB function
        matlab_function = getattr(matlab_engine, function_name)
        
        # Execute the function
        if nargout == 0:
            matlab_function(*args, nargout=0, **kwargs)
            result = None
        else:
            result = matlab_function(*args, nargout=nargout, **kwargs)
        
        # Convert MATLAB arrays to Python lists if needed
        if result is not None:
            try:
                if hasattr(result, 'tolist'):
                    result = result.tolist()
                elif isinstance(result, tuple):
                    result = [r.tolist() if hasattr(r, 'tolist') else r for r in result]
            except:
                pass  # Keep result as-is if conversion fails

        print(
            f"[MATLAB HTTP] Execution completed: function={function_name}, result_summary={str(result)[:200] if result is not None else 'None'}",
            flush=True
        )
        
        return jsonify({
            'status': 'success',
            'result': result
        }), 200
        
    except matlab.engine.MatlabExecutionError as e:
        print(f"[MATLAB HTTP] MATLAB execution error for function {function_name}: {e}", flush=True)
        return jsonify({
            'status': 'error',
            'message': f'MATLAB execution error: {str(e)}'
        }), 500
    except Exception as e:
        print(f"[MATLAB HTTP] Unexpected error for function {locals().get('function_name', 'unknown')}: {e}", flush=True)
        print(traceback.format_exc(), flush=True)
        return jsonify({
            'status': 'error',
            'message': f'Unexpected error: {str(e)}',
            'traceback': traceback.format_exc()
        }), 500

@app.route('/eval', methods=['POST'])
def eval_code():
    """
    Evaluate MATLAB code
    
    Request body:
    {
        "code": "matlab code string",
        "nargout": 1  // optional, number of output arguments (default 0)
    }
    
    Response:
    {
        "status": "success" | "error",
        "result": <result_value>,  // if status is success and nargout > 0
        "message": <error_message>  // if status is error
    }
    """
    try:
        data = request.get_json()
        
        if not data or 'code' not in data:
            return jsonify({
                'status': 'error',
                'message': 'Missing required field: code'
            }), 400
        
        code = data['code']
        nargout = data.get('nargout', 0)

        print(
            f"[MATLAB HTTP] Received eval request: nargout={nargout}, code_snippet={code[:100]}...",
            flush=True
        )
        
        # Execute the code
        if nargout == 0:
            matlab_engine.eval(code, nargout=0)
            result = None
        else:
            result = matlab_engine.eval(code, nargout=nargout)
        
        # Convert MATLAB arrays to Python lists if needed
        if result is not None:
            try:
                if hasattr(result, 'tolist'):
                    result = result.tolist()
            except:
                pass

        print("[MATLAB HTTP] Eval completed", flush=True)
        
        return jsonify({
            'status': 'success',
            'result': result
        }), 200
        
    except matlab.engine.MatlabExecutionError as e:
        print(f"[MATLAB HTTP] MATLAB execution error during eval: {e}", flush=True)
        return jsonify({
            'status': 'error',
            'message': f'MATLAB execution error: {str(e)}'
        }), 500
    except Exception as e:
        print(f"[MATLAB HTTP] Unexpected error during eval: {e}", flush=True)
        print(traceback.format_exc(), flush=True)
        return jsonify({
            'status': 'error',
            'message': f'Unexpected error: {str(e)}',
            'traceback': traceback.format_exc()
        }), 500

def shutdown_matlab():
    """Ensure MATLAB engine is closed cleanly on exit."""
    global matlab_engine
    if matlab_engine is not None:
        try:
            print("Shutting down MATLAB engine...", flush=True)
            matlab_engine.quit()
        except Exception:
            pass

# Register shutdown handler
atexit.register(shutdown_matlab)

if __name__ == "__main__":
    # Initialize MATLAB before starting the Flask server
    try:
        initialize_matlab()
    except Exception as e:
        print(f"Failed to initialize MATLAB: {e}", file=sys.stderr, flush=True)
        sys.exit(1)
    
    # Start Flask server
    port = int(os.getenv('MATLAB_HTTP_PORT', '9090'))
    print(f"Starting MATLAB HTTP API server on port {port}...", flush=True)
    app.run(host='0.0.0.0', port=port, debug=False, threaded=True)
