"""
Solution 4 — Nested Model with Validation — Vaccine Profile
=============================================================
Builds nested Pydantic models (DosePlan, ManufacturerInfo,
VaccineProfile) with a field validator, and uses structured
output to parse LLM-generated vaccine information.

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

# Step 1: Define the nested models
class DosePlan(BaseModel):
    route: str = Field(description="Administration route, e.g., intramuscular")
    dose_value: float = Field(description="Dose amount")
    dose_unit: str = Field(description="Dose unit, e.g., mL, mg")

class ManufacturerInfo(BaseModel):
    name: str = Field(description="Manufacturer company name")
    country: Optional[str] = Field(None, description="Country of origin")

class VaccineProfile(BaseModel):
    name: str = Field(description="Vaccine name")
    indications: List[str] = Field(description="What the vaccine protects against")
    schedule: DosePlan = Field(description="Dosing information")
    manufacturer: ManufacturerInfo = Field(description="Manufacturer details")

    # Step 2: Validator — dose_value must be positive
    @field_validator("schedule")
    def dose_must_be_positive(cls, v):
        if v.dose_value <= 0:
            raise ValueError("dose_value must be > 0")
        return v

# Step 3: Generate a paragraph about a vaccine
response = llm.invoke(
    "Write a short paragraph about the Hepatitis B vaccine including: "
    "name, indications, administration route, dose with unit, "
    "and manufacturer details."
)
vaccine_text = response.content
print("Generated text:")
print(vaccine_text)
print()

# Step 4: Parse with structured output
structured_llm = llm.with_structured_output(VaccineProfile)
result = structured_llm.invoke(f"Extract vaccine profile from this text: {vaccine_text}")

# Step 5: Print the validated object
print("Parsed VaccineProfile:")
print(result)
print()
print(f"Name: {result.name}")
print(f"Indications: {result.indications}")
print(f"Schedule: {result.schedule.dose_value} {result.schedule.dose_unit} via {result.schedule.route}")
print(f"Manufacturer: {result.manufacturer.name} ({result.manufacturer.country})")
