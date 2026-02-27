"""
Solution 3 — Union Classification for Route-Specific Therapies
===============================================================
Defines TopicalTherapy and InhaledTherapy models with a
TherapyWrapper that classifies which type of therapy a text
describes using optional fields and a model_validator.

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

# Step 1: Define the therapy models
class TopicalTherapy(BaseModel):
    name: str = Field(description="Name of the topical therapy")
    formulation: str = Field(description="e.g., cream, ointment, gel")
    site_of_application: str = Field(description="Where it is applied")

class InhaledTherapy(BaseModel):
    name: str = Field(description="Name of the inhaled therapy")
    device_type: str = Field(description="e.g., HFA inhaler, DPI")
    dose_mcg_per_puff: Optional[float] = Field(None, description="Dose per puff in mcg")

# Step 2: Create the wrapper model with a model_validator
class TherapyWrapper(BaseModel):
    """Holds EITHER a topical OR an inhaled therapy, never both."""
    topical: Optional[TopicalTherapy] = Field(None, description="Fill ONLY if the therapy is topical")
    inhaled: Optional[InhaledTherapy] = Field(None, description="Fill ONLY if the therapy is inhaled")

    @model_validator(mode="after")
    def exactly_one_set(self):
        if self.topical is not None and self.inhaled is not None:
            raise ValueError("Only one of topical or inhaled should be set, not both")
        if self.topical is None and self.inhaled is None:
            raise ValueError("At least one of topical or inhaled must be set")
        return self

    def get_therapy_type(self) -> str:
        if self.topical is not None:
            return "topical"
        return "inhaled"

# Step 3: Parse two inputs
structured_llm = llm.with_structured_output(TherapyWrapper)

topical_text = "Miconazole 2% cream applied twice daily to affected skin area."
inhaled_text = "Albuterol HFA inhaler 100 mcg per puff used for wheezing."

prompt = "Extract the therapy info. Fill ONLY the topical field if topical, or ONLY the inhaled field if inhaled: {text}"

topical_result = structured_llm.invoke(prompt.format(text=topical_text))
inhaled_result = structured_llm.invoke(prompt.format(text=inhaled_text))

# Step 4: Print results
print(f"Topical input → detected type: {topical_result.get_therapy_type()}")
print(f"  Details: {topical_result}")
print()
print(f"Inhaled input → detected type: {inhaled_result.get_therapy_type()}")
print(f"  Details: {inhaled_result}")
