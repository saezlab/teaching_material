"""
Exercise 2 — Local vs Cloud Comparison
========================================
Compare a local Ollama model with a cloud-based Gemini model.

Your tasks:
  1. Create a local LLM (Ollama) and a cloud LLM (Gemini)
  2. Define a shared question
  3. Invoke both models, measure response time
  4. Print and compare the results

Usage (from the playground/ root):
    uv run lecture_2/4_ollama/exercise_2_local_vs_cloud.py
"""

import time
from dotenv import load_dotenv
load_dotenv()

from langchain.chat_models import init_chat_model

# --- Configuration ---
# MODEL_NAME = "gemini-3-flash-preview"
# MODEL_PROVIDER = "google_genai"
MODEL_NAME = "my-model"
MODEL_PROVIDER = "openai"
BASE_URL = "https://myung-courageous-gino.ngrok-free.dev/v1"
API_KEY = "not-needed"

# ── Step 1: Create the local model ──────────────────────────────
# TODO: Use init_chat_model() with:
#   - model="qwen3:0.6b"
#   - model_provider="ollama"
#   - temperature=0.2

local_llm = None  # <-- replace this

# ── Step 2: Create the cloud model ──────────────────────────────
# TODO: Use init_chat_model() with:
#   - model=MODEL_NAME
#   - model_provider=MODEL_PROVIDER
#   - base_url=BASE_URL
#   - api_key=API_KEY
#   - temperature=0.2

cloud_llm = None  # <-- replace this

# ── Step 3: Define the question ──────────────────────────────────
question = "Explain how vaccines work in simple terms."

# ── Step 4: Query both models and compare ────────────────────────
# TODO: For each model (local_llm and cloud_llm):
#   1. Record the start time with time.time()
#   2. Invoke the model with the question
#   3. Record the end time
#   4. Print the model name, response text (response.text), and elapsed time
#   5. Print a separator line between models
#
# Hint: You can put them in a list of tuples:
#   models = [("Local (qwen3:0.6b)", local_llm), ("Cloud (Gemini)", cloud_llm)]
#   for name, llm in models:
#       ...

# ── Feedback ────────────────────────────────────────────────────
if local_llm is None or cloud_llm is None:
    print("Exercise not yet completed — edit this file to fill in the TODOs above.")
