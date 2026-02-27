"""
Exercise 1 — LangChain Tools
==============================
Define biology-themed tools with @tool, bind them to an LLM,
execute tool calls, and return results for a final answer.

Your tasks:
  1. Define a melting_temperature tool (Wallace rule)
  2. Define a reverse_complement tool
  3. Bind both tools to the LLM
  4. Execute tool calls from the model's response
  5. Return results to the model for a final narrative answer

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
# TODO: Define a function called melting_temperature that takes a
# sequence (str) and returns a float.
#
# Use the Wallace rule for short oligonucleotides:
#   Tm = 2 * (number of A + T bases) + 4 * (number of G + C bases)
#
# Convert the sequence to uppercase first.
# Don't forget the @tool decorator and a clear docstring!


# ── Step 2: Define the reverse_complement tool ───────────
# TODO: Define a function called reverse_complement that takes a
# sequence (str) and returns a str.
#
# The complement mapping is: A<->T, C<->G
# The reverse complement is the complement read backwards.
# For example: ATGC -> complement is TACG -> reverse is GCAT
#
# Convert the sequence to uppercase first.
# Don't forget the @tool decorator and a clear docstring!


# ── Step 3: Bind tools to the LLM ───────────────────────
# TODO: Create llm_with_tools by calling llm.bind_tools([...])
# with your two tools.

llm_with_tools = None  # Replace this

if llm_with_tools is None:
    print("\nExercise not yet completed — fill in the TODOs above!")
    exit()

# ── Invoke the model ────────────────────────────────────

prompt = ChatPromptTemplate.from_messages([
    ("system", "You are a helpful molecular biology assistant. Use tools when needed."),
    ("human", "{question}")
])

chain = prompt | llm_with_tools

question = "What is the melting temperature of the DNA sequence ATGCGATCG?"
print(f"Question: {question}\n")

ai_msg = chain.invoke({"question": question})

# Check if the model made any tool calls
tool_calls = getattr(ai_msg, "tool_calls", [])
print(f"Number of tool calls: {len(tool_calls)}")

if not tool_calls:
    print("No tool calls were made. The model answered directly:")
    print(ai_msg.content)
else:
    # ── Step 4: Execute tool calls ───────────────────────
    # TODO: Create a dictionary mapping tool names to tool objects.
    # Then loop through tool_calls, look up each tool by name,
    # and invoke it with the provided arguments.
    #
    # Example:
    #   tools_by_name = {"melting_temperature": melting_temperature, ...}
    #   for call in tool_calls:
    #       result = tools_by_name[call["name"]].invoke(call["args"])
    #       print(f"Tool: {call['name']}, Result: {result}")

    pass  # Replace this with your code

    # ── Step 5: Return results to the model ──────────────
    # TODO: Build a message list and send tool results back to the
    # model so it can write a final answer.
    #
    # The message list should contain:
    #   1. The original prompt messages
    #   2. The AI message with tool calls (ai_msg)
    #   3. A ToolMessage for each tool result
    #
    # Then invoke the base llm (not llm_with_tools) with these messages.
    #
    # Note: The guide used HumanMessage(...) to build the list manually.
    # Here we use prompt.format_messages() which does the same thing
    # but includes the system message automatically.
    #
    # Example:
    #   messages = prompt.format_messages(question=question)
    #   messages.append(ai_msg)
    #   for call in tool_calls:
    #       result = tools_by_name[call["name"]].invoke(call["args"])
    #       messages.append(ToolMessage(content=str(result), tool_call_id=call["id"]))
    #   final = llm.invoke(messages)
    #   print(final.content)

    pass  # Replace this with your code
