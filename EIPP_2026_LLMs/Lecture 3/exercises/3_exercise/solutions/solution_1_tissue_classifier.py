"""
Solution 1 — LLM as a Tissue Classifier
=========================================
Uses a simple Pydantic model with a Literal field and
with_structured_output() to classify histological tissue
excerpts into tissue types.

Usage (from the playground/ root):
    uv run lecture_3/3_structured_outputs/exercise_1_tissue_classifier.py
"""

from dotenv import load_dotenv
load_dotenv()

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pydantic")

from langchain.chat_models import init_chat_model
from pydantic import BaseModel, Field
from typing import List, Literal

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

# Step 1: Define the Pydantic model with a Literal for allowed tissue origins
class TissueClassification(BaseModel):
    organ: Literal["heart", "liver", "brain", "bone", "skin"] = Field(
        description="The organ this tissue sample comes from"
    )
    key_features: List[str] = Field(
        description="Key histological features found in the excerpt"
    )

# Step 2: The five excerpts to classify
excerpts = [
    "The biopsy revealed striated cells with intercalated discs and central nuclei, "
    "consistent with myocardial tissue. Branching fibers were observed throughout the sample.",

    "Histological examination showed hepatocytes arranged in lobular architecture "
    "with central veins and portal triads. Bile canaliculi and Kupffer cells were identified "
    "within the sinusoidal lining.",

    "The sample contained neurons with extensive dendritic arborization and myelinated axon "
    "bundles. Glial cells, including astrocytes and oligodendrocytes, were abundant.",

    "Cross-section analysis revealed osteons with concentric lamellae surrounding Haversian "
    "canals. Osteocytes were observed within lacunae connected by canaliculi.",

    "The specimen displayed stratified squamous epithelium with a distinct keratinized layer. "
    "Underlying dermal papillae and collagen bundles were visible in the supporting connective tissue.",
]

# Step 3: Classify each excerpt
classifier = llm.with_structured_output(TissueClassification)

print("Tissue Classification Results")
print("=" * 50)

for i, excerpt in enumerate(excerpts, 1):
    result = classifier.invoke(
        f"Classify which organ this tissue sample comes from: {excerpt}"
    )
    print(f"\nExcerpt {i}: {excerpt[:60]}...")
    print(f"  Organ:         {result.organ}")
    print(f"  Key features:  {result.key_features}")
