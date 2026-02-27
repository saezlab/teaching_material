"""
Solution — Exercise 1: First OpenRouter Call
"""

import os
from dotenv import load_dotenv
load_dotenv()

from langchain.chat_models import init_chat_model

# Step 1: Create the model
llm = init_chat_model(
    model="google/gemma-3-27b-it:free",
    model_provider="openai",
    base_url="https://openrouter.ai/api/v1",
    api_key=os.environ["OPENROUTER_API_KEY"],
)

# Step 2: Invoke the model
response = llm.invoke("What are the main differences between DNA and RNA?")

# Step 3: Print the response and token usage
print("--- Response ---")
print(response.text)
print("\n--- Token usage ---")
print(response.usage_metadata)
