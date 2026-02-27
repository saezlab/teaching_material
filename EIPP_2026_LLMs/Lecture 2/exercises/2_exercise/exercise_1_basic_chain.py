"""
Exercise 1 — Basic Chain Construction
======================================
Build a LangChain chain that explains scientific concepts simply.

Your tasks:
  1. Create a ChatPromptTemplate with a system message and a {concept} placeholder
  2. Build a chain: prompt | llm | StrOutputParser()
  3. Invoke the chain with {"concept": "photosynthesis"}
  4. Print the result

Usage (from the playground/ root):
    uv run lecture_2/2_langchain/exercise_1_basic_chain.py
"""

from dotenv import load_dotenv
load_dotenv()

from langchain.chat_models import init_chat_model
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser

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

# ── Step 1: Create the prompt template ──────────────────────────────
# TODO: Use ChatPromptTemplate.from_messages() with:
#   - A ("system", "...") message telling the model to act as a friendly
#     science teacher who explains complex topics in simple terms
#   - A ("human", "{concept}") message

prompt = None  # <-- replace this

# ── Step 2: Build the chain ─────────────────────────────────────────
# TODO: Connect prompt | llm | StrOutputParser()

chain = None  # <-- replace this

# ── Step 3: Invoke and print ────────────────────────────────────────
# TODO: Call chain.invoke({"concept": "photosynthesis"}) and print the result

# ── Feedback ────────────────────────────────────────────────────────
if chain is None:
    print("Exercise not yet completed — edit this file to fill in the TODOs above.")
