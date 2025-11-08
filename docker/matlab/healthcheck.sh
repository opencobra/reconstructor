#!/bin/bash
# Health check script for MATLAB container
# This script verifies that MATLAB Engine API is running and accessible

set -e

echo "Checking MATLAB Engine API health..."

# Check if Python and MATLAB Engine are available
python3 -c "import matlab.engine" 2>/dev/null || {
    echo "ERROR: MATLAB Engine for Python not available"
    exit 1
}

# Check if shared MATLAB sessions exist
SESSION_COUNT=$(python3 -c "import matlab.engine; print(len(matlab.engine.find_matlab()))" 2>/dev/null || echo "0")

if [ "$SESSION_COUNT" -eq "0" ]; then
    echo "ERROR: No shared MATLAB sessions found"
    exit 1
fi

echo "SUCCESS: Found $SESSION_COUNT MATLAB session(s)"

# List available sessions
python3 -c "import matlab.engine; print('Available sessions:', matlab.engine.find_matlab())"

exit 0
