"""
Exercise 1 — ReAct Agent with Biology Tools
=============================================
Build a ReAct agent using LangGraph that uses biology-themed tools
to answer questions about DNA sequences.

Your tasks:
  1. Define the AgentState TypedDict
  2. Implement the agent_node function
  3. Implement the route_after_agent function
  4. Build and compile the ReAct graph

Usage (from the playground/ root):
    uv run lecture_5/1_react_agent/exercise_react_agent.py
"""

from typing import Annotated
from typing_extensions import TypedDict
from langchain_core.messages import HumanMessage, SystemMessage
from langchain_core.tools import tool
from langgraph.graph import StateGraph, END
from langgraph.graph.message import add_messages
from langgraph.prebuilt import ToolNode
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


# ── Tools (provided — do not edit) ──────────────────────────────

@tool
def gc_content(sequence: str) -> float:
    """Calculate the GC content of a DNA sequence as a percentage."""
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    return round((gc_count / len(sequence)) * 100, 2)


@tool
def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    sequence = sequence.upper()
    return "".join(complement[base] for base in reversed(sequence))


TOOLS = [gc_content, reverse_complement]


# ── TODO 1: Define the AgentState ───────────────────────────────
# Create a TypedDict with a "messages" key.
# Use Annotated[list, add_messages] as the type.
# ─────────────────────────────────────────────────────────────────


# ── TODO 2: Define the agent_node ───────────────────────────────
# Write a function agent_node(state: AgentState) that:
#   1. Creates an LLM with init_chat_model(
#          model=MODEL_NAME,
#          model_provider=MODEL_PROVIDER,
#          base_url=BASE_URL,
#          api_key=API_KEY,
#      )
#   2. Binds the tools: llm_with_tools = llm.bind_tools(TOOLS)
#   3. Invokes the LLM: ai_msg = llm_with_tools.invoke(state["messages"])
#   4. Returns {"messages": [ai_msg]}
# ─────────────────────────────────────────────────────────────────


# ── TODO 3: Define the routing function ─────────────────────────
# Write a function route_after_agent(state: AgentState) -> str:
#   - Get the last message: last = state["messages"][-1]
#   - Check for tool calls: getattr(last, "tool_calls", [])
#   - Return "tools" if there are tool calls, "end" otherwise
# ─────────────────────────────────────────────────────────────────


# ── TODO 4: Build the ReAct graph ───────────────────────────────
# 1. Create a StateGraph(AgentState)
# 2. Add node "agent" with agent_node
# 3. Add node "tools" with ToolNode(TOOLS)
# 4. Set "agent" as entry point
# 5. Add conditional edges from "agent":
#        workflow.add_conditional_edges(
#            "agent", route_after_agent, {"tools": "tools", "end": END}
#        )
# 6. Add edge from "tools" back to "agent" (this creates the loop!)
# 7. Compile: graph = workflow.compile()
#
# graph = ...
# ─────────────────────────────────────────────────────────────────


graph = None  # Replace this with your compiled graph

# ── Run the agent (provided — do not edit) ──────────────────────
if __name__ == "__main__":
    if graph is None:
        print("\nExercise not yet completed — fill in the TODOs above!")
        exit()

    system_prompt = (
        "You are a helpful biology assistant. Use your tools when "
        "they are relevant to answer questions about DNA sequences."
    )

    questions = [
        "What is the GC content of ATGCGATCG?",
        "What is the reverse complement of ATGCGATCG?",
        "What is the GC content of the reverse complement of ATGCGATCG?",
    ]

    for question in questions:
        print(f"\nQuestion: {question}")
        result = graph.invoke({
            "messages": [
                SystemMessage(content=system_prompt),
                HumanMessage(content=question),
            ]
        })
        print(f"Answer: {result['messages'][-1].content}")
        print("-" * 60)
