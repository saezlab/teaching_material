"""
Exercise 3 — Open Targets MCP
===============================
Connect to the remote Open Targets MCP server and use an agent
to answer a biology question about a gene of your choice.

Your tasks:
  1. Connect to the Open Targets MCP server
  2. Load and inspect the available tools
  3. Create an agent with create_agent
  4. Ask a biology question and print the answer

Usage (from the playground/ root):
    uv run lecture_4/3_open_targets/exercise_open_targets.py
"""

from dotenv import load_dotenv
load_dotenv()

import asyncio

from mcp import ClientSession
from mcp.client.streamable_http import streamablehttp_client
from langchain_mcp_adapters.tools import load_mcp_tools
from langchain.agents import create_agent
from langchain.chat_models import init_chat_model

# --- Configuration ---
# MODEL_NAME = "gemini-3-flash-preview"
# MODEL_PROVIDER = "google_genai"
MODEL_NAME = "my-model"
MODEL_PROVIDER = "openai"
BASE_URL = "https://myung-courageous-gino.ngrok-free.dev/v1"
API_KEY = "not-needed"

async def main():
    server_url = "https://mcp.platform.opentargets.org/mcp"

    # ── Step 1: Connect to the MCP server ────────────────
    # TODO: Use streamablehttp_client to connect to server_url,
    # then open a ClientSession and initialize it.
    #
    # The pattern is:
    #   async with streamablehttp_client(server_url) as (read, write, _):  # The _ is a third value we don't need
    #       async with ClientSession(read, write) as session:
    #           await session.initialize()
    #           ... (steps 2-4 go here, inside this block)

    # ── Step 2: Load and inspect tools ───────────────────
    # TODO: Load the MCP tools and print their names:
    #   tools = await load_mcp_tools(session)
    #   print(f"Found {len(tools)} tools:")
    #   for t in tools:
    #       print(f"  - {t.name}")

    # ── Step 3: Create the agent ─────────────────────────
    # TODO: Create the LLM and agent:
    #   llm = init_chat_model(
    #       model=MODEL_NAME,
    #       model_provider=MODEL_PROVIDER,
    #       base_url=BASE_URL,
    #       api_key=API_KEY,
    #       temperature=0.2,
    #   )
    #   agent = create_agent(model=llm, tools=tools)

    # ── Step 4: Ask a question ───────────────────────────
    # TODO: Ask a biology question. Pick a gene that interests you!
    #
    # Examples:
    #   "What diseases are associated with TP53?"
    #   "What drugs target EGFR?"
    #   "What is known about the BRCA1 gene?"
    #
    #   result = await agent.ainvoke({
    #       "messages": [("user", "YOUR QUESTION HERE")]
    #   })
    #   print(result["messages"][-1].content)

    print("Exercise not yet completed — fill in the TODOs above!")


if __name__ == "__main__":
    asyncio.run(main())
