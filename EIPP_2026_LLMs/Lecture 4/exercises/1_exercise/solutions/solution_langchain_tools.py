"""
Solution 1 — LangChain Tools
==============================
Defines biology-themed tools, binds them to an LLM, executes tool calls,
and returns results for a final narrative answer.

Usage (from the playground/ root):
    uv run lecture_4/1_langchain_tools/exercise_langchain_tools.py
"""

from dotenv import load_dotenv
load_dotenv()

from langchain.chat_models import init_chat_model
from langchain_core.tools import tool
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.messages import ToolMessage

# --- Configuration ---
# MODEL_NAME = "gemini-3-flash-preview"
# MODEL_PROVIDER = "google_genai"
MODEL_NAME = "my-model"
MODEL_PROVIDER = "openai"
BASE_URL = "https://myung-courageous-gino.ngrok-free.dev/v1"
API_KEY = "not-needed"

# ── Initialize the model ─────────────────────────────────

llm = init_chat_model(
    model=MODEL_NAME,
    model_provider=MODEL_PROVIDER,
    base_url=BASE_URL,
    api_key=API_KEY,
    temperature=0.2,
)

# ── Step 1: Define the melting_temperature tool ──────────

@tool
def melting_temperature(sequence: str) -> float:
    """Estimate the melting temperature of a short DNA sequence using the Wallace rule.

    The Wallace rule: Tm = 2*(A+T count) + 4*(G+C count).
    Best for oligonucleotides of 14 bases or fewer.

    Args:
        sequence: a DNA sequence string (e.g. "ATGCGATCG")

    Returns:
        the estimated melting temperature in degrees Celsius
    """
    seq = sequence.upper()
    at_count = seq.count("A") + seq.count("T")
    gc_count = seq.count("G") + seq.count("C")
    return 2.0 * at_count + 4.0 * gc_count

# ── Step 2: Define the reverse_complement tool ───────────

@tool
def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Uses standard base pairing: A<->T, C<->G.

    Args:
        sequence: a DNA sequence string (e.g. "ATGCGATCG")

    Returns:
        the reverse complement as a string
    """
    complement_map = {"A": "T", "T": "A", "C": "G", "G": "C"}
    seq = sequence.upper()
    complement = "".join(complement_map.get(base, base) for base in seq)
    return complement[::-1]

# ── Step 3: Bind tools to the LLM ───────────────────────

llm_with_tools = llm.bind_tools([melting_temperature, reverse_complement])

# ── Step 4: Invoke the model ────────────────────────────

prompt = ChatPromptTemplate.from_messages([
    ("system", "You are a helpful molecular biology assistant. Use tools when needed."),
    ("human", "{question}")
])

chain = prompt | llm_with_tools

question = "What is the melting temperature of the DNA sequence ATGCGATCG?"
print(f"Question: {question}\n")

ai_msg = chain.invoke({"question": question})

tool_calls = getattr(ai_msg, "tool_calls", [])
print(f"Number of tool calls: {len(tool_calls)}")

if not tool_calls:
    print("No tool calls were made. The model answered directly:")
    print(ai_msg.content)
else:
    # Execute tool calls
    tools_by_name = {
        "melting_temperature": melting_temperature,
        "reverse_complement": reverse_complement,
    }

    tool_results = []
    for call in tool_calls:
        name = call["name"]
        args = call["args"]
        result = tools_by_name[name].invoke(args)
        tool_results.append({"name": name, "result": result, "id": call["id"]})
        print(f"Tool: {name}, Args: {args}, Result: {result}")

    # Return results to the model for a final answer
    messages = prompt.format_messages(question=question)
    messages.append(ai_msg)

    for tr in tool_results:
        messages.append(
            ToolMessage(
                content=str(tr["result"]),
                name=tr["name"],
                tool_call_id=tr["id"],
            )
        )

    final_answer = llm.invoke(messages)
    print(f"\nFinal answer:\n{final_answer.content}")
