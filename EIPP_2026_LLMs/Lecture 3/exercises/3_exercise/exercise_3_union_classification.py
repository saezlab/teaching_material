"""
Exercise 3 — Union Classification for Route-Specific Therapies
===============================================================
Define two Pydantic models for different therapy routes and use
a wrapper with optional fields so the LLM correctly classifies
the therapy type.

Your tasks:
  1. Define TopicalTherapy and InhaledTherapy models
  2. Create a TherapyWrapper with optional fields and a model_validator
  3. Parse two inputs (topical and inhaled)
  4. Print the detected type for each

Usage (from the playground/ root):
    uv run lecture_3/3_structured_outputs/exercise_3_union_classification.py
"""

from dotenv import load_dotenv
load_dotenv()

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pydantic")

from langchain.chat_models import init_chat_model
from pydantic import BaseModel, Field, model_validator
from typing import Optional

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

# ── Step 1: Define the therapy models ────────────────────
# TODO: Define a class TopicalTherapy(BaseModel) with:
#   - name: str               — name of the topical therapy
#   - formulation: str         — e.g., cream, ointment, gel
#   - site_of_application: str — where it is applied

# TODO: Define a class InhaledTherapy(BaseModel) with:
#   - name: str                       — name of the inhaled therapy
#   - device_type: str                — e.g., HFA inhaler, DPI
#   - dose_mcg_per_puff: Optional[float] — dose per puff in mcg (default None)


# ── Step 2: Create the wrapper model with a model_validator ────────────────────
# TODO: Define a class TherapyWrapper(BaseModel) with:
#   - topical: Optional[TopicalTherapy] = Field(None, description="Fill ONLY if the therapy is topical")
#   - inhaled: Optional[InhaledTherapy] = Field(None, description="Fill ONLY if the therapy is inhaled")
#
# Add a @model_validator(mode="after") that checks:
#   - If both topical and inhaled are set → raise ValueError
#   - If neither is set → raise ValueError
#
# Also add a helper method get_therapy_type() that returns "topical" or "inhaled"
# based on which field is not None.


# ── Step 3: Parse two inputs ─────────────────────────────
# TODO: Create a structured LLM using llm.with_structured_output(TherapyWrapper).
# Then parse these two texts:

topical_text = "Miconazole 2% cream applied twice daily to affected skin area."
inhaled_text = "Albuterol HFA inhaler 100 mcg per puff used for wheezing."

# Use a clear prompt that tells the LLM to fill only one field:
# prompt = "Extract the therapy info. Fill ONLY the topical field if topical, or ONLY the inhaled field if inhaled: {text}"
# Store results in topical_result and inhaled_result.

topical_result = None
inhaled_result = None


# ── Step 4: Print results ────────────────────────────────
# TODO: Print the detected type for each result.
# For example:
#   print(f"Topical input → detected type: {topical_result.get_therapy_type()}")
#   print(f"  Details: {topical_result}")
#   print()
#   print(f"Inhaled input → detected type: {inhaled_result.get_therapy_type()}")
#   print(f"  Details: {inhaled_result}")


# ── Feedback ─────────────────────────────────────────────
if topical_result is None or inhaled_result is None:
    print("Exercise not yet completed — fill in the TODOs above!")
