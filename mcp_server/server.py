"""Compatibility shim: the server lives in src/constructor_mcp/server.py.

Kept so existing MCP client configs pointing at `python server.py` keep
working. New installs should prefer the packaged entry point:

    uvx constructor-mcp          # or: pip install constructor-mcp
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent / 'src'))

from constructor_mcp.server import main  # noqa: E402

if __name__ == '__main__':
    main()
