#!/usr/bin/env bash
# Build the one-click Claude Desktop bundle: mcp_server/build/constructor-mcp.mcpb
#
# Requires: python3 with pip (override with PIP=...), and the mcpb CLI
# (npm install -g @anthropic-ai/mcpb).
set -euo pipefail
cd "$(dirname "$0")"

PIP="${PIP:-python3 -m pip}"

rm -rf build/mcpb
mkdir -p build/mcpb/server

cp mcpb/manifest.json build/mcpb/manifest.json
cp src/constructor_mcp/server.py build/mcpb/server/main.py

# Vendor the runtime dependencies so users need nothing but a system python3.
$PIP install --quiet --target build/mcpb/server/lib 'mcp>=2,<3' 'httpx>=0.27'

mcpb validate build/mcpb/manifest.json
mcpb pack build/mcpb build/constructor-mcp.mcpb

echo
echo "Bundle ready: $(pwd)/build/constructor-mcp.mcpb"
echo "Install: open it with Claude Desktop (double-click), fill in URL + API key."
