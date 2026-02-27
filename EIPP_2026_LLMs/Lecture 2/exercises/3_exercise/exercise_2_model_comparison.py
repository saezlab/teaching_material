"""
Exercise 2 — Model Comparison
==============================
Compare two free OpenRouter models on the same question.

Your tasks:
  1. Create a helper function to build OpenRouter models
  2. Call multiple models with the same question
  3. Measure and compare response time and quality

Usage (from the playground/ root):
    uv run lecture_2/3_openrouter/exercise_2_model_comparison.py
"""

import os
import time
from dotenv import load_dotenv
load_dotenv()

from langchain.chat_models import init_chat_model

# ── Step 1: Create a helper function ─────────────────────────────
# TODO: Write a function make_openrouter_llm(model_name) that returns
#       an init_chat_model instance configured for OpenRouter.
#
# Hint:
#   def make_openrouter_llm(model_name):
#       return init_chat_model(
#           model=model_name,
#           model_provider="openai",
#           base_url="https://openrouter.ai/api/v1",
#           api_key=os.environ["OPENROUTER_API_KEY"],
#       )

def make_openrouter_llm(model_name):
    pass  # <-- replace this

# ── Step 2: Define models and question ───────────────────────────
# TODO: Create a list of 2 free model names to compare.
# Browse https://openrouter.ai/models?pricing=free to find free models.
# Example: ["google/gemma-3-27b-it:free", "deepseek/deepseek-r1-0528:free"]

models = []  # <-- replace this

question = "Explain how CRISPR-Cas9 works in simple terms."

# ── Step 3: Compare the models ───────────────────────────────────
# TODO: For each model in the list:
#   1. Create the LLM with make_openrouter_llm()
#   2. Record the start time with time.time()
#   3. Invoke the model with the question
#   4. Record the end time
#   5. Print:
#      - Model name
#      - Response text (response.text)
#      - Token usage (response.usage_metadata)
#      - Elapsed time (end - start)
#   6. Print a separator line between models

# ── Feedback ────────────────────────────────────────────────────
if not models:
    print("Exercise not yet completed — edit this file to fill in the TODOs above.")
