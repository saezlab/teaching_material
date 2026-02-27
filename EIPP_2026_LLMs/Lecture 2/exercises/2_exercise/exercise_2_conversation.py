"""
Exercise 2 — Message Types and Conversation
============================================
Simulate a doctor-patient conversation about medication side effects.

Your tasks:
  1. Create a SystemMessage for the doctor persona
  2. Build a conversation with at least 2 exchanges
  3. Print each response
  4. Print the type of each message at the end

Usage (from the playground/ root):
    uv run lecture_2/2_langchain/exercise_2_conversation.py
"""

from dotenv import load_dotenv
load_dotenv()

from langchain.chat_models import init_chat_model
from langchain_core.messages import SystemMessage, HumanMessage, AIMessage

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

# ── Step 1: Create the system message ──────────────────────────────
# TODO: Create a SystemMessage that makes the model act as a
#       knowledgeable but empathetic doctor

system_msg = None  # <-- replace this

# ── Step 2: First exchange ─────────────────────────────────────────
# TODO:
#   - Create a HumanMessage asking about side effects of a medication
#   - Build a messages list: [system_msg, human_msg]
#   - Invoke the model with the messages list
#   - Store the response (it is already an AIMessage)
#   - Print the response with response.text
#
# Hint — the pattern for each exchange is:
#   messages.append(human_msg)
#   response = llm.invoke(messages)
#   messages.append(response)

# ── Step 3: Second exchange (follow-up) ────────────────────────────
# TODO:
#   - Create another HumanMessage with a follow-up question
#   - Append both the previous AI response and this new question to the list
#   - Invoke the model with the full conversation history
#   - Print the response with response.text

# ── Step 4: Print message types ────────────────────────────────────
# TODO: Loop through the conversation history and print the type of each message
# Example output:
#   <class 'langchain_core.messages.system.SystemMessage'>
#   <class 'langchain_core.messages.human.HumanMessage'>
#   ...

# ── Feedback ────────────────────────────────────────────────────────
if system_msg is None:
    print("Exercise not yet completed — edit this file to fill in the TODOs above.")
