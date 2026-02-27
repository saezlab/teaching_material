"""
Exercise 3 — Multi-Specialist Consultation System
==================================================
Create specialist chains using partial() and MessagesPlaceholder.

Your tasks:
  1. Create a base prompt with partial() for different specialists
  2. Use MessagesPlaceholder for conversation history
  3. Build a make_specialist_chain() function
  4. Show different specialists answering the same question
  5. Print token usage after each response

Usage (from the playground/ root):
    uv run lecture_2/2_langchain/exercise_3_multi_specialist.py
"""

from dotenv import load_dotenv
load_dotenv()

from langchain.chat_models import init_chat_model
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain_core.output_parsers import StrOutputParser
from langchain_core.messages import HumanMessage, AIMessage

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

# ── Step 1: Create the base prompt template ─────────────────────────
# TODO: Create a ChatPromptTemplate with:
#   - A system message containing a {specialty} placeholder
#     (e.g., "You are a specialist in {specialty}. Answer patient questions
#      from the perspective of your specialty.")
#   - A MessagesPlaceholder("history") for conversation history
#   - A ("human", "{question}") for the current question

base_prompt = None  # <-- replace this

# ── Step 2: Create specialist prompts using partial() ───────────────
# TODO: Use base_prompt.partial(specialty="...") to create prompts for
#       2 specialties (e.g., cardiology and neurology)

# cardiology_prompt = base_prompt.partial(specialty="cardiology")
# neurology_prompt = ...

# ── Step 3: Create a make_specialist_chain() function ───────────────
# TODO: Write a function that takes a specialty name and returns a chain.
#       Hint: use base_prompt.partial(specialty=name) | llm
#       Note: do NOT add StrOutputParser() here — we need the full AIMessage
#       to access usage_metadata for token tracking.

def make_specialist_chain(specialty: str):
    pass  # <-- replace this

# ── Step 4: Ask the same question to two specialists ────────────────
# TODO:
#   - Pick a patient question (e.g., "I have been experiencing chest pain
#     and shortness of breath. What could be causing this?")
#   - Invoke each specialist chain with {"history": [], "question": your_question}
#   - Save the responses (e.g., in a dict) so you can reuse one for step 5
#   - Print each specialist's response with response.text
#   - Print the token usage from response.usage_metadata

# ── Step 5: Demonstrate conversation history ────────────────────────
# TODO:
#   - Pick one specialist (e.g., cardiology) and reuse its response from step 4
#   - Build a history list from the step 4 exchange:
#       history = [HumanMessage(content=patient_question), responses["cardiology"]]
#   - Ask a follow-up question (e.g., "What tests would you recommend?")
#     passing the history so the specialist remembers context
#   - Print the follow-up response

# ── Feedback ────────────────────────────────────────────────────────
if base_prompt is None:
    print("Exercise not yet completed — edit this file to fill in the TODOs above.")
