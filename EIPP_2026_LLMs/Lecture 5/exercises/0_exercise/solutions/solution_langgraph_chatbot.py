"""
Solution — Exercise 0: LangGraph Chatbot
==========================================
A simple chatbot using LangGraph's StateGraph.

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


# ── Step 1: Define the State ────────────────────────────────────

class State(TypedDict):
    messages: Annotated[list, add_messages]


# ── Step 2: Define the call_llm node ────────────────────────────

def call_llm(state: State):
    llm = init_chat_model(
        model=MODEL_NAME,
        model_provider=MODEL_PROVIDER,
        base_url=BASE_URL,
        api_key=API_KEY,
        temperature=0.2,
    )
    ai_msg = llm.invoke(state["messages"])
    return {"messages": [ai_msg]}


# ── Step 3: Build and compile the graph ─────────────────────────

workflow = StateGraph(State)
workflow.add_node("llm", call_llm)
workflow.set_entry_point("llm")
workflow.add_edge("llm", END)
graph = workflow.compile()


# ── Chat loop ───────────────────────────────────────────────────
if __name__ == "__main__":
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
