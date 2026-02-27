"""
Exercise 2 — Structured Output for a Dietary Supplement
========================================================
Define a Pydantic model for a dietary supplement and use
LangChain's with_structured_output() to extract structured
data from LLM-generated text.

Your tasks:
  1. Define a Supplement Pydantic model with the fields listed below
  2. Generate a short paragraph about "vitamin D3"
  3. Use with_structured_output(Supplement) to parse the text
  4. Print and verify the parsed object

Usage (from the playground/ root):
    uv run lecture_3/3_structured_outputs/exercise_2_supplement.py
"""

from dotenv import load_dotenv
load_dotenv()

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pydantic")

from langchain.chat_models import init_chat_model
from pydantic import BaseModel, Field
from typing import List, Optional

# --- Configuration ---
# MODEL_NAME = "gemini-3-flash-preview"
# MODEL_PROVIDER = "google_genai"
MODEL_NAME = "my-model"
MODEL_PROVIDER = "openai"
BASE_URL = "https://myung-courageous-gino.ngrok-free.dev/v1"
API_KEY = "not-needed"

# Initialize the model
llm = init_chat_model(
    model=MODEL_NAME,
    model_provider=MODEL_PROVIDER,
    base_url=BASE_URL,
    api_key=API_KEY,
    temperature=0.2,
)

# ── Step 1: Define the Supplement model ──────────────────
# TODO: Define a class called Supplement that inherits from BaseModel.
# It should have these fields (use Field() with a description for each):
#   - common_names: List[str]           — common names for the supplement
#   - active_components: List[str]      — active chemical components
#   - uses: List[str]                   — therapeutic or health uses
#   - typical_adult_dose_mg: Optional[float] — typical adult dose in mg (default None)
#   - forms: List[str]                  — available forms (tablet, capsule, liquid, etc.)


# ── Step 2: Generate a paragraph about vitamin D3 ───────
# TODO: Use llm.invoke() to ask the LLM to write a short paragraph about
# vitamin D3 covering: common names, active components, uses, typical
# adult dose in mg, and available forms.
# Store the response text in a variable called supplement_text.

supplement_text = None


# ── Step 3: Parse with structured output ─────────────────
# TODO: Create a structured LLM using llm.with_structured_output(Supplement).
# Then invoke it with a prompt that asks it to extract supplement info
# from supplement_text.
# Store the result in a variable called result.

result = None


# ── Step 4: Print and verify ─────────────────────────────
# TODO: Print the result and its individual fields to verify they look correct.
# For example:
#   print(f"Common names: {result.common_names}")
#   print(f"Active components: {result.active_components}")
#   ... and so on for each field.


# ── Feedback ─────────────────────────────────────────────
if result is None:
    print("Exercise not yet completed — fill in the TODOs above!")
