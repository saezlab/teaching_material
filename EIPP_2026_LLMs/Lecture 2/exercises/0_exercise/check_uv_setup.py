"""
UV Setup Verification Script
=============================
Run this script inside your uv-managed project to verify
that your environment is correctly set up.

Usage:
    uv run check_uv_setup.py
"""

import sys
import os

REQUIRED_PACKAGES = ["cowsay"]
EXPECTED_PYTHON = (3, 11)


def check_virtualenv():
    """Check that we are running inside a virtual environment."""
    in_venv = sys.prefix != sys.base_prefix
    venv_path = os.environ.get("VIRTUAL_ENV", "")
    if in_venv:
        print(f"  Virtual environment: {sys.prefix}")
    else:
        print("  WARNING: Not running inside a virtual environment.")
        print("  Make sure you activated it or use 'uv run' to run this script.")
    return in_venv


def check_python_version():
    """Check that the Python version matches the expected one."""
    major, minor = sys.version_info[:2]
    match = (major, minor) == EXPECTED_PYTHON
    print(f"  Python version: {major}.{minor} ", end="")
    if match:
        print(f"(expected {EXPECTED_PYTHON[0]}.{EXPECTED_PYTHON[1]})")
    else:
        print(f"(expected {EXPECTED_PYTHON[0]}.{EXPECTED_PYTHON[1]}, MISMATCH)")
    return match


def check_packages():
    """Check that all required packages can be imported."""
    all_ok = True
    for pkg in REQUIRED_PACKAGES:
        try:
            mod = __import__(pkg)
            version = getattr(mod, "__version__", "unknown")
            print(f"  {pkg}: {version}")
        except ImportError:
            print(f"  {pkg}: NOT FOUND")
            all_ok = False
    return all_ok


def main():
    print("=" * 50)
    print("UV Setup Verification")
    print("=" * 50)

    results = {}

    print("\n[1/3] Checking virtual environment...")
    results["virtualenv"] = check_virtualenv()

    print("\n[2/3] Checking Python version...")
    results["python"] = check_python_version()

    print("\n[3/3] Checking required packages...")
    results["packages"] = check_packages()

    print("\n" + "=" * 50)
    print("Results:")
    print("=" * 50)

    all_passed = True
    for name, passed in results.items():
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False

    print()
    if all_passed:
        import cowsay
        cowsay.cow("Moo! All checks passed! Your UV environment is set up correctly.")
    else:
        print("Some checks failed. Review the output above and fix the issues.")

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
