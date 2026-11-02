"""Smoke test for the Constructor MCP server (US4).

Spawns server.py over stdio with the MCP client SDK, checks the tool surface,
and — if CONSTRUCTOR_API_KEY is set and the API is reachable — performs one
live resolve_metabolite round trip.

Run:  python test_smoke.py
"""
import asyncio
import os
import sys

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

EXPECTED_TOOLS = {
    'resolve_metabolite', 'create_reaction', 'save_reaction',
    'list_reactions', 'get_reaction', 'fetch_vmh_reaction',
    'fetch_rhea_reaction',
}

SERVER = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'server.py')


async def main():
    params = StdioServerParameters(
        command=sys.executable, args=[SERVER],
        env={**os.environ},
    )
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            tools = (await session.list_tools()).tools
            names = {tool.name for tool in tools}
            assert names == EXPECTED_TOOLS, f'tool surface mismatch: {names}'
            for tool in tools:
                text = (tool.description or '').lower()
                assert 'submit' not in tool.name.lower()
                assert text, f'tool {tool.name} has no description'
            print(f'OK: {len(tools)} tools, surface matches contract')

            if os.environ.get('CONSTRUCTOR_API_KEY'):
                result = await session.call_tool(
                    'resolve_metabolite',
                    {'identifier': 'O', 'id_type': 'smiles'})
                text = result.content[0].text
                if text.startswith('Error [connection]'):
                    print('SKIP live call: API not reachable '
                          f'({os.environ.get("CONSTRUCTOR_BASE_URL", "http://localhost:8000")})')
                else:
                    assert '"formula": "H2O"' in text, text
                    print('OK: live resolve_metabolite returned H2O')
            else:
                print('SKIP live call: CONSTRUCTOR_API_KEY not set')


if __name__ == '__main__':
    asyncio.run(main())
