"""
Exercise 2 — MCP Client
=========================
Connect to the MCP math server, load its tools, and use an
agent to answer math questions.

Your tasks:
  1. Set up the server connection with StdioServerParameters
  2. Load MCP tools with load_mcp_tools
  3. Create an agent with create_agent
  4. Test with a question that uses your new tools

Usage (from the playground/ root):
    uv run lecture_4/2_mcp/exercise_mcp_client.py
"""

from dotenv import load_dotenv
load_dotenv()

import asyncio
import os

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
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
    # ── Step 1: Set up the server connection ─────────────
    # __file__ is a special variable containing this script's path. We use it to find math_server.py in the same folder.
    server_path = os.path.join(os.path.dirname(__file__), "math_server.py")

    # TODO: Create StdioServerParameters to connect to the math server.
    # Use "uv" as the command and ["run", server_path] as the args,
    # so the server runs inside the same virtual environment.
    #
    # Example:
    #   server_params = StdioServerParameters(
    #       command="uv",
    #       args=["run", server_path],
    #   )

    server_params = None  # Replace this

    # ── Step 2: Connect and load tools ───────────────────
    # TODO: Use stdio_client and ClientSession to connect to the server,
    # initialize the session, and load tools.
    #
    # The pattern is:
    #   async with stdio_client(server_params) as (read, write):
    #       async with ClientSession(read, write) as session:
    #           await session.initialize()
    #           tools = await load_mcp_tools(session)
    #           ... (steps 3 and 4 go here, inside this block)

    # ── Step 3: Create the agent ─────────────────────────
    # TODO: Create an LLM and an agent:
    #   llm = init_chat_model(
    #       model=MODEL_NAME,
    #       model_provider=MODEL_PROVIDER,
    #       base_url=BASE_URL,
    #       api_key=API_KEY,
    #       temperature=0.2,
    #   )
    #   agent = create_agent(model=llm, tools=tools)

    # ── Step 4: Ask a question ───────────────────────────
    # TODO: Invoke the agent with a question that uses your new tools.
    #
    # Example:
    #   result = await agent.ainvoke({
    #       "messages": [("user", "What is 2 to the power of 10?")]
    #   })
    #   print(result["messages"][-1].content)

    print("Exercise not yet completed — fill in the TODOs above!")


if __name__ == "__main__":
    asyncio.run(main())
