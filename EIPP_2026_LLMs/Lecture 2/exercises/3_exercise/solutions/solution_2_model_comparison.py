"""
Solution — Exercise 2: Model Comparison
"""

import os
import time
from dotenv import load_dotenv
load_dotenv()

from langchain.chat_models import init_chat_model


def make_openrouter_llm(model_name):
    """Create a LangChain chat model that uses OpenRouter."""
    return init_chat_model(
        model=model_name,
        model_provider="openai",
        base_url="https://openrouter.ai/api/v1",
        api_key=os.environ["OPENROUTER_API_KEY"],
    )


# Define models to compare (free models on OpenRouter)
models = [
    "google/gemma-3-27b-it:free",
    "deepseek/deepseek-r1-0528:free",
]

question = "Explain how CRISPR-Cas9 works in simple terms."

# Compare each model
for model_name in models:
    print(f"=== {model_name} ===")

    llm = make_openrouter_llm(model_name)

    start = time.time()
    response = llm.invoke(question)
    elapsed = time.time() - start

    print(response.text)
    print(f"\nTokens: {response.usage_metadata}")
    print(f"Time:   {elapsed:.1f}s")
    print()
