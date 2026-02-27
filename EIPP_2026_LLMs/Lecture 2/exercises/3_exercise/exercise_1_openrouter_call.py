"""
Exercise 1 — First OpenRouter Call
===================================
Connect to OpenRouter via LangChain and call a free model.

Your tasks:
  1. Create an LLM with init_chat_model pointing to OpenRouter
  2. Invoke the model with a biology question
  3. Print the response and token usage

Usage (from the playground/ root):
    uv run lecture_2/3_openrouter/exercise_1_openrouter_call.py
"""

import os
from dotenv import load_dotenv
load_dotenv()

from langchain.chat_models import init_chat_model

# ── Step 1: Create the model ─────────────────────────────────────
# TODO: Create an LLM with init_chat_model():
#   - model="google/gemma-3-27b-it:free"
#   - model_provider="openai"
#   - base_url="https://openrouter.ai/api/v1"
#   - api_key=os.environ["OPENROUTER_API_KEY"]

llm = None  # <-- replace this

# ── Step 2: Invoke the model ─────────────────────────────────────
# TODO: Call llm.invoke() with a biology question
#   e.g., "What are the main differences between DNA and RNA?"
# Store the result in a variable called `response`

# ── Step 3: Print the response and token usage ───────────────────
# TODO: Print the response text with response.text
# TODO: Print the token usage with response.usage_metadata

# ── Feedback ────────────────────────────────────────────────────
if llm is None:
    print("Exercise not yet completed — edit this file to fill in the TODOs above.")
