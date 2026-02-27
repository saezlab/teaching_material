"""
Exercise 0 — LangGraph Chatbot
================================
Build a simple chatbot using LangGraph.
Follow the TODOs below to complete the exercise.

Your tasks:
  1. Define the State TypedDict
  2. Implement the call_llm node function
  3. Build and compile the graph

Usage (from the playground/ root):
    uv run lecture_5/0_langgraph_chatbot/exercise_langgraph_chatbot.py
"""

from typing import Annotated
from typing_extensions import TypedDict
from langchain_core.messages import HumanMessage
from langgraph.graph import StateGraph, START, END
from langgraph.graph.message import add_messages
from langchain.chat_models import init_chat_model
from dotenv import load_dotenv

load_dotenv()

# --- Configuration ---
# MODEL_NAME = "gemini-3-flash-preview"
# MODEL_PROVIDER = "google_genai"
MODEL_NAME = "my-model"
MODEL_PROVIDER = "openai"
BASE_URL = "https://myung-courageous-gino.ngrok-free.dev/v1"
API_KEY = "not-needed"


# ── TODO 1: Define the State ────────────────────────────────────
# Create a TypedDict called State with a single key "messages".
# The type should be Annotated[list, add_messages] so that
# LangGraph appends new messages instead of overwriting the list.
#
# Example:
#   class MyState(TypedDict):
#       some_key: Annotated[list, add_messages]
# ─────────────────────────────────────────────────────────────────


# ── TODO 2: Define the call_llm node ────────────────────────────
# Write a function call_llm(state: State) that:
#   1. Creates an LLM with init_chat_model(
#          model=MODEL_NAME,
#          model_provider=MODEL_PROVIDER,
#          base_url=BASE_URL,
#          api_key=API_KEY,
#      )
#   2. Invokes the LLM with state["messages"]
#   3. Returns {"messages": [ai_msg]}
# ─────────────────────────────────────────────────────────────────


# ── TODO 3: Build and compile the graph ─────────────────────────
# 1. Create a StateGraph with your State class
# 2. Add a node called "llm" with call_llm as the function
# 3. Set "llm" as the entry point
# 4. Add an edge from "llm" to END
# 5. Compile the graph and store it in the variable `graph`
#
# graph = ...
# ─────────────────────────────────────────────────────────────────


graph = None  # Replace this with your compiled graph

# ── Chat loop (provided — do not edit) ──────────────────────────
if __name__ == "__main__":
    if graph is None:
        print("\nExercise not yet completed — fill in the TODOs above!")
        exit()

    print("LangGraph Chatbot (type 'quit' to exit)")
    print("-" * 40)

    state = None
    while True:
        user_input = input("\nYou: ")
        if user_input.lower() == "quit":
            break

        if state is None:
            state = graph.invoke(
                {"messages": [HumanMessage(content=user_input)]}
            )
        else:
            state["messages"].append(HumanMessage(content=user_input))
            state = graph.invoke(state)

        print(f"\nAssistant: {state['messages'][-1].content}")
