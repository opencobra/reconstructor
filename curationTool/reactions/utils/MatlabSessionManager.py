import json
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
        # Get the base directory (two levels up from this script)
        base_dir = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                '..',
                '..',
                '..'))
        config_path = os.path.join(base_dir, 'config.json')

        if not os.path.exists(config_path):
            raise FileNotFoundError(
                f"Configuration file {config_path} not found.")

        with open(config_path, 'r') as config_file:
            config = json.load(config_file)

        script_directories = config.get('script_directories', [])
        cobra_path = config.get('cobra_path', '')

        for script_directory in script_directories:
            self.engine.addpath(script_directory, nargout=0)
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
