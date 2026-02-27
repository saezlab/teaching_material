"""
Solution — Exercise 1: ReAct Agent with Biology Tools
======================================================
A ReAct agent using LangGraph that uses biology-themed tools
to answer questions about DNA sequences.

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


# ── Step 1: Define the AgentState ───────────────────────────────

class AgentState(TypedDict):
    messages: Annotated[list, add_messages]


# ── Step 2: Define the agent_node ───────────────────────────────

def agent_node(state: AgentState):
    llm = init_chat_model(
        model=MODEL_NAME,
        model_provider=MODEL_PROVIDER,
        base_url=BASE_URL,
        api_key=API_KEY,
        temperature=0.2,
    )
    llm_with_tools = llm.bind_tools(TOOLS)
    ai_msg = llm_with_tools.invoke(state["messages"])
    return {"messages": [ai_msg]}


# ── Step 3: Define the routing function ─────────────────────────

def route_after_agent(state: AgentState) -> str:
    last = state["messages"][-1]
    if getattr(last, "tool_calls", []):
        return "tools"
    return "end"


# ── Step 4: Build the ReAct graph ───────────────────────────────

workflow = StateGraph(AgentState)
workflow.add_node("agent", agent_node)
workflow.add_node("tools", ToolNode(TOOLS))
workflow.set_entry_point("agent")
workflow.add_conditional_edges(
    "agent", route_after_agent, {"tools": "tools", "end": END}
)
workflow.add_edge("tools", "agent")
graph = workflow.compile()


# ── Run the agent ───────────────────────────────────────────────
if __name__ == "__main__":
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
