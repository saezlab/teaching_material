"""
Gemini API Key Verification Script
====================================
Verifies that your GOOGLE_API_KEY environment variable is set
and working by making a simple request to the Gemini API.

Uses only Python standard library (no external dependencies).

Usage:
    uv run check_api_key.py
"""

import json
import os
import sys
import urllib.request
import urllib.error

API_URL = "https://generativelanguage.googleapis.com/v1beta/models/gemini-3-flash-preview:generateContent"
TEST_PROMPT = "Reply with exactly: Hello, the API key is working!"


def main():
    print("=" * 50)
    print("Gemini API Key Verification")
    print("=" * 50)

    # --- Check that the environment variable is set ---
    print("\n[1/2] Checking GOOGLE_API_KEY environment variable...")
    api_key = os.environ.get("GOOGLE_API_KEY")

    if not api_key:
        print("  FAIL: GOOGLE_API_KEY is not set.")
        print()
        print("  Set it in your terminal before running this script:")
        print()
        print("    export GOOGLE_API_KEY='your-api-key-here'")
        print()
        sys.exit(1)

    print(f"  Found: {api_key[:8]}...{api_key[-4:]}")

    # --- Test the key with a real API call ---
    print("\n[2/2] Testing API call to Gemini...")

    payload = json.dumps({
        "contents": [{"parts": [{"text": TEST_PROMPT}]}],
    }).encode("utf-8")

    request = urllib.request.Request(
        f"{API_URL}?key={api_key}",
        data=payload,
        headers={"Content-Type": "application/json"},
        method="POST",
    )

    try:
        with urllib.request.urlopen(request) as response:
            result = json.loads(response.read().decode("utf-8"))
            text = result["candidates"][0]["content"]["parts"][0]["text"]
            print(f"  Gemini says: {text.strip()}")
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", errors="replace")
        print(f"  FAIL: HTTP {e.code}")
        print(f"  {body[:300]}")
        sys.exit(1)
    except Exception as e:
        print(f"  FAIL: {e}")
        sys.exit(1)

    print("\n" + "=" * 50)
    print("All checks passed! Your API key is working.")
    print("=" * 50)


if __name__ == "__main__":
    main()
