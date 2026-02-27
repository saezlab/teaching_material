"""
Ollama Verification Script
===========================
Checks that Ollama is installed, running, and has the required model.

Usage (from the playground/ root):
    uv run lecture_2/4_ollama/check_ollama.py
"""

import subprocess
import sys


def check_ollama_installed():
    """Check that the ollama command is available."""
    try:
        result = subprocess.run(
            ["ollama", "--version"],
            capture_output=True, text=True, timeout=10,
        )
        if result.returncode == 0:
            print(f"[OK] Ollama is installed: {result.stdout.strip()}")
            return True
        else:
            print("[FAIL] Ollama command returned an error.")
            return False
    except FileNotFoundError:
        print("[FAIL] Ollama is not installed. Visit https://ollama.com/download")
        return False
    except subprocess.TimeoutExpired:
        print("[FAIL] Ollama command timed out.")
        return False


def check_ollama_running():
    """Check that the Ollama service is reachable."""
    try:
        import urllib.request
        req = urllib.request.Request("http://localhost:11434/api/tags")
        with urllib.request.urlopen(req, timeout=5) as resp:
            if resp.status == 200:
                print("[OK] Ollama service is running.")
                return True
    except Exception:
        pass
    print("[FAIL] Ollama service is not running.")
    print("       On macOS: start the Ollama app.")
    print("       On Linux: run 'ollama serve' in a separate terminal.")
    return False


def check_model_available(model_name="qwen3:0.6b"):
    """Check that the required model is pulled."""
    try:
        result = subprocess.run(
            ["ollama", "list"],
            capture_output=True, text=True, timeout=10,
        )
        if model_name.split(":")[0] in result.stdout:
            print(f"[OK] Model '{model_name}' is available.")
            return True
        else:
            print(f"[FAIL] Model '{model_name}' not found.")
            print(f"       Run: ollama pull {model_name}")
            return False
    except Exception:
        print("[FAIL] Could not check available models.")
        return False


if __name__ == "__main__":
    print("Checking Ollama setup...\n")
    all_ok = True
    all_ok = check_ollama_installed() and all_ok
    all_ok = check_ollama_running() and all_ok
    all_ok = check_model_available() and all_ok

    print()
    if all_ok:
        print("All checks passed! You are ready for the Ollama exercises.")
    else:
        print("Some checks failed. Please fix the issues above before continuing.")
        sys.exit(1)
