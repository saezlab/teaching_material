"""
Exercise 4 — Nested Model with Validation — Vaccine Profile
=============================================================
Build nested Pydantic models with a field validator for a
vaccine profile.

Your tasks:
  1. Define DosePlan, ManufacturerInfo, and VaccineProfile models
  2. Add a @field_validator ensuring dose_value > 0
  3. Generate a paragraph about a vaccine
  4. Parse it with with_structured_output(VaccineProfile)
  5. Print the validated object

Usage (from the playground/ root):
    uv run lecture_3/3_structured_outputs/exercise_4_nested_vaccine.py
"""

from dotenv import load_dotenv
load_dotenv()

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pydantic")

from langchain.chat_models import init_chat_model
from pydantic import BaseModel, Field, field_validator
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

# ── Step 1: Define the nested models ────────────────────
# TODO: Define a class DosePlan(BaseModel) with:
#   - route: str       — administration route (e.g., intramuscular)
#   - dose_value: float — the dose amount
#   - dose_unit: str    — the unit (e.g., mL, mg)

# TODO: Define a class ManufacturerInfo(BaseModel) with:
#   - name: str             — manufacturer company name
#   - country: Optional[str] — country of origin (default None)

# TODO: Define a class VaccineProfile(BaseModel) with:
#   - name: str                    — vaccine name
#   - indications: List[str]       — what the vaccine protects against
#   - schedule: DosePlan           — dosing information (nested model)
#   - manufacturer: ManufacturerInfo — manufacturer details (nested model)


# ── Step 2: Add a validator ──────────────────────────────
# TODO: Inside VaccineProfile, add a @field_validator("schedule")
# that checks if schedule.dose_value > 0. If not, raise a ValueError.
# Example:
#
#   @field_validator("schedule")
#   def dose_must_be_positive(cls, v):
#       if v.dose_value <= 0:
#           raise ValueError("dose_value must be > 0")
#       return v


# ── Step 3: Generate a paragraph about a vaccine ────────
# TODO: Use llm.invoke() to ask the LLM to write a short paragraph about
# the Hepatitis B vaccine, including: name, indications, administration
# route, dose with unit, and manufacturer details.
# Store the response text in a variable called vaccine_text.

vaccine_text = None


# ── Step 4: Parse with structured output ─────────────────
# TODO: Create a structured LLM using llm.with_structured_output(VaccineProfile).
# Then invoke it to extract the vaccine profile from vaccine_text.
# Store the result in a variable called result.

result = None


# ── Step 5: Print the validated object ───────────────────
# TODO: Print the result and its nested fields. For example:
#   print(f"Name: {result.name}")
#   print(f"Indications: {result.indications}")
#   print(f"Schedule: {result.schedule.dose_value} {result.schedule.dose_unit} via {result.schedule.route}")
#   print(f"Manufacturer: {result.manufacturer.name} ({result.manufacturer.country})")


# ── Feedback ─────────────────────────────────────────────
if result is None:
    print("Exercise not yet completed — fill in the TODOs above!")
