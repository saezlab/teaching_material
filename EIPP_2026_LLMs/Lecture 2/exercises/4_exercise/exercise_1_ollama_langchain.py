"""
Exercise 1 — Ollama with LangChain
====================================
Use LangChain to query a locally-deployed Ollama model.

Your tasks:
  1. Create an LLM with init_chat_model (model_provider="ollama")
  2. Create a ChatPromptTemplate with a system message and {question} placeholder
  3. Build a chain: prompt | llm | StrOutputParser()
  4. Invoke the chain with a biology question
  5. Print the result

Usage (from the playground/ root):
    uv run lecture_2/4_ollama/exercise_1_ollama_langchain.py
"""

from langchain.chat_models import init_chat_model
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser

# ── Step 1: Create the model ─────────────────────────────────────
# TODO: Use init_chat_model() with:
#   - model="qwen3:0.6b"
#   - model_provider="ollama"
#   - temperature=0.2

llm = None  # <-- replace this

# ── Step 2: Create the prompt template ───────────────────────────
# TODO: Use ChatPromptTemplate.from_messages() with:
#   - A ("system", "...") message telling the model to act as a helpful
#     biology tutor who explains topics clearly
#   - A ("human", "{question}") message

prompt = None  # <-- replace this

# ── Step 3: Build the chain ──────────────────────────────────────
# TODO: Connect prompt | llm | StrOutputParser()

chain = None  # <-- replace this

# ── Step 4: Invoke and print ─────────────────────────────────────
# TODO: Call chain.invoke({"question": "What is photosynthesis?"})
#       and print the result

# ── Feedback ────────────────────────────────────────────────────
if chain is None:
    print("Exercise not yet completed — edit this file to fill in the TODOs above.")
