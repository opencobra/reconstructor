import os
import matlab.engine


class MatlabSessionManager:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(MatlabSessionManager, cls).__new__(cls)
            try:
                cls._instance.engine = matlab.engine.start_matlab()
                cls._instance._setup_cobra_toolbox()
            except Exception as e:
                print(f"Failed to start MATLAB session: {e}")
                cls._instance = None
        return cls._instance

    def _setup_cobra_toolbox(self):
        # Read configuration from environment variables
        script_directories_env = os.getenv('SCRIPT_DIRECTORIES', '')
        script_directories = [d.strip() for d in script_directories_env.split(',') if d.strip()]
        cobra_path = os.getenv('COBRA_PATH', '')

        for script_directory in script_directories:
            self.engine.addpath(script_directory, nargout=0)
        if cobra_path:
            self.engine.addpath(cobra_path, nargout=0)
            self.engine.eval("initCobraToolbox(0)", nargout=0)

    def execute(self, command, *args, **kwargs):
        try:
            if hasattr(self.engine, command):
                matlab_function = getattr(self.engine, command)
                result = matlab_function(*args, **kwargs)
                return {'status': 'success', 'result': result}
            else:
                return {
                    'status': 'error',
                    'message': f'Command {command} not found'}
        except matlab.engine.MatlabExecutionError as e:
            return {'status': 'error', 'message': str(e)}

    def quit(self):
        if self.engine:
            self.engine.quit()
            self.__class__._instance = None
