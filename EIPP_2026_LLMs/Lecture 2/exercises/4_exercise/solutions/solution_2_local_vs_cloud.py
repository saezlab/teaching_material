"""
Solution — Exercise 2: Local vs Cloud Comparison
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

# Step 1: Create both models
local_llm = init_chat_model(
    model="qwen3:0.6b",
    model_provider="ollama",
    temperature=0.2,
)

cloud_llm = init_chat_model(
    model=MODEL_NAME,
    model_provider=MODEL_PROVIDER,
    base_url=BASE_URL,
    api_key=API_KEY,
    temperature=0.2,
)

# Step 2: Define the question
question = "Explain how vaccines work in simple terms."

# Step 3: Query both models and compare
models = [
    ("Local (qwen3:0.6b)", local_llm),
    ("Cloud (EBI endpoint)", cloud_llm),
]

for name, llm in models:
    print(f"=== {name} ===")

    start = time.time()
    response = llm.invoke(question)
    elapsed = time.time() - start

    print(response.text)
    print(f"\nTime: {elapsed:.1f}s")
    print()
