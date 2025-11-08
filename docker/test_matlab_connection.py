#!/usr/bin/env python3
"""
Test script for MATLAB Session Manager

This script tests the connection to the MATLAB container and verifies
that basic MATLAB operations work correctly.

Usage:
    python test_matlab_connection.py
"""

import os
import sys

# Add the project to the path
sys.path.insert(0, '/app/curationTool')

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'reactions_project.settings')

import django
django.setup()

def test_matlab_connection():
    """Test the MATLAB connection and basic operations."""
    
    print("=" * 60)
    print("MATLAB Session Manager Connection Test")
    print("=" * 60)
    
    # Check environment
    print(f"\nEnvironment Configuration:")
    print(f"  MATLAB_REMOTE_ENABLED: {os.getenv('MATLAB_REMOTE_ENABLED', 'not set')}")
    print(f"  MATLAB_SESSION_NAME: {os.getenv('MATLAB_SESSION_NAME', 'not set')}")
    
    # Import the session manager
    try:
        if os.getenv('MATLAB_REMOTE_ENABLED', 'false').lower() == 'true':
            from reactions.utils.MatlabSessionManagerRemote import MatlabSessionManager
            print("  Using: Remote MATLAB Session Manager")
        else:
            from reactions.utils.MatlabSessionManager import MatlabSessionManager
            print("  Using: Local MATLAB Session Manager")
    except ImportError as e:
        print(f"\n❌ Failed to import MatlabSessionManager: {e}")
        return False
    
    # Initialize the session manager
    print("\n" + "-" * 60)
    print("Initializing MATLAB Session...")
    print("-" * 60)
    
    try:
        matlab = MatlabSessionManager()
        print("✓ Session manager initialized successfully")
    except Exception as e:
        print(f"❌ Failed to initialize session manager: {e}")
        return False
    
    # Check connection status
    if hasattr(matlab, 'is_connected'):
        if matlab.is_connected:
            print("✓ MATLAB engine is connected and responsive")
        else:
            print("❌ MATLAB engine is not responsive")
            return False
    
    # Test basic operations
    print("\n" + "-" * 60)
    print("Testing Basic MATLAB Operations...")
    print("-" * 60)
    
    tests = [
        {
            'name': 'Simple arithmetic',
            'command': 'eval',
            'args': ('2 + 2',),
            'kwargs': {'nargout': 1},
            'expected': 4.0
        },
        {
            'name': 'Square root',
            'command': 'sqrt',
            'args': (16.0,),
            'kwargs': {'nargout': 1},
            'expected': 4.0
        },
        {
            'name': 'Array creation',
            'command': 'ones',
            'args': (2, 2),
            'kwargs': {'nargout': 1},
            'expected': None  # Just check it doesn't error
        },
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            result = matlab.execute(
                test['command'],
                *test['args'],
                **test['kwargs']
            )
            
            if result['status'] == 'success':
                if test['expected'] is None or result['result'] == test['expected']:
                    print(f"✓ {test['name']}: PASSED")
                    passed += 1
                else:
                    print(f"❌ {test['name']}: FAILED (unexpected result)")
                    print(f"   Expected: {test['expected']}")
                    print(f"   Got: {result['result']}")
                    failed += 1
            else:
                print(f"❌ {test['name']}: FAILED")
                print(f"   Error: {result.get('message', 'Unknown error')}")
                failed += 1
        except Exception as e:
            print(f"❌ {test['name']}: FAILED (exception)")
            print(f"   Error: {e}")
            failed += 1
    
    # Test custom function (if generateVMHMetAbbr exists)
    print("\n" + "-" * 60)
    print("Testing Custom Functions...")
    print("-" * 60)
    
    try:
        result = matlab.execute('generateVMHMetAbbr', 'glucose')
        if result['status'] == 'success':
            print(f"✓ generateVMHMetAbbr: PASSED")
            print(f"   Result: {result['result']}")
            passed += 1
        else:
            print(f"⚠ generateVMHMetAbbr: Not available or failed")
            print(f"   Message: {result.get('message', 'Unknown')}")
    except Exception as e:
        print(f"⚠ generateVMHMetAbbr: Not available")
        print(f"   Error: {e}")
    
    # Summary
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    
    if failed == 0:
        print("\n✓ All tests passed! MATLAB connection is working correctly.")
        return True
    else:
        print(f"\n❌ {failed} test(s) failed. Please check the configuration.")
        return False

if __name__ == '__main__':
    success = test_matlab_connection()
    sys.exit(0 if success else 1)
