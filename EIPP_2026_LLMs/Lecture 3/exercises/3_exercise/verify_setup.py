"""
Lecture 3 — Setup Verification
================================
Run this script to check that your environment is ready for the
structured outputs exercises.

Usage (from the playground/ root):
    uv run lecture_3/3_structured_outputs/verify_setup.py
"""

import sys

errors = []

# 1. Check Python version
if sys.version_info < (3, 11):
    errors.append(f"Python 3.11+ required, but you have {sys.version}")

# 2. Check pydantic
try:
    import pydantic
    print(f"[OK] pydantic {pydantic.__version__}")
except ImportError:
    errors.append("pydantic not installed — run: uv add pydantic")

# 3. Check langchain
try:
    from langchain.chat_models import init_chat_model
    print("[OK] langchain")
except ImportError:
    errors.append("langchain not installed — run: uv add langchain")

# 4. Check python-dotenv
try:
    from dotenv import load_dotenv
    print("[OK] python-dotenv")
except ImportError:
    errors.append("python-dotenv not installed — run: uv add python-dotenv")

# 5. Check LLM endpoint
try:
    from langchain.chat_models import init_chat_model
    llm = init_chat_model(
        model="my-model",
        model_provider="openai",
        base_url="https://myung-courageous-gino.ngrok-free.dev/v1",
        api_key="not-needed",
        temperature=0.2,
    )
    response = llm.invoke("Say OK")
    if response and response.content:
        print("[OK] LLM endpoint is reachable")
    else:
        errors.append("LLM endpoint returned an empty response")
except Exception as e:
    errors.append(
        f"LLM endpoint unreachable ({e}). "
        "If the EBI endpoint is down, uncomment the Gemini lines in the exercise scripts instead."
    )

# Summary
print()
if errors:
    print("PROBLEMS FOUND:")
    for err in errors:
        print(f"  - {err}")
    sys.exit(1)
else:
    print("All checks passed — you're ready for the Lecture 3 exercises!")
