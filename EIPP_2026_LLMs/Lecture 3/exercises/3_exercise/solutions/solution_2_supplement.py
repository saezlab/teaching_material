"""
Solution 2 — Structured Output for a Dietary Supplement
========================================================
Defines a Supplement Pydantic model and uses LangChain's
with_structured_output() to extract structured data from
LLM-generated text about vitamin D3.

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

# Step 1: Define the Supplement model
class Supplement(BaseModel):
    common_names: List[str] = Field(description="Common names for the supplement")
    active_components: List[str] = Field(description="Active chemical components")
    uses: List[str] = Field(description="Therapeutic or health uses")
    typical_adult_dose_mg: Optional[float] = Field(None, description="Typical adult dose in mg")
    forms: List[str] = Field(description="Available forms (e.g., tablet, capsule, liquid)")

# Step 2: Generate a paragraph about vitamin D3
response = llm.invoke(
    "Write a short paragraph about vitamin D3 covering: "
    "common names, active components, uses, typical adult dose in mg, "
    "and available forms."
)
supplement_text = response.content
print("Generated text:")
print(supplement_text)
print()

# Step 3: Parse with structured output
structured_llm = llm.with_structured_output(Supplement)
result = structured_llm.invoke(f"Extract supplement info from this text: {supplement_text}")

# Step 4: Print and verify
print("Parsed Supplement object:")
print(result)
print()
print(f"Common names: {result.common_names}")
print(f"Active components: {result.active_components}")
print(f"Uses: {result.uses}")
print(f"Typical dose (mg): {result.typical_adult_dose_mg}")
print(f"Forms: {result.forms}")
