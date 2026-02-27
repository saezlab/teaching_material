"""
Solution 2 — MCP Client
=========================
Connects to the extended MCP math server and uses an agent
to answer math questions.

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
    # Path to the math server (same directory as this script)
    server_path = os.path.join(os.path.dirname(__file__), "math_server.py")

    # Set up the server connection using stdio transport
    server_params = StdioServerParameters(
        command="uv",
        args=["run", server_path],
    )

    async with stdio_client(server_params) as (read, write):
        async with ClientSession(read, write) as session:
            # Initialize the MCP connection and discover tools
            await session.initialize()
            tools = await load_mcp_tools(session)

            print(f"Discovered {len(tools)} tools:")
            for t in tools:
                print(f"  - {t.name}")
            print()

            # Create the LLM and agent
            llm = init_chat_model(
                model=MODEL_NAME,
                model_provider=MODEL_PROVIDER,
                base_url=BASE_URL,
                api_key=API_KEY,
                temperature=0.2,
            )
            agent = create_agent(model=llm, tools=tools)

            # Test questions
            questions = [
                "What is 2 to the power of 10?",
                "What is 100 minus 37?",
                "What is 3 + 5, and then raise that result to the power of 2?",
            ]

            for question in questions:
                print(f"Question: {question}")
                result = await agent.ainvoke({
                    "messages": [("user", question)]
                })
                print(f"Answer: {result['messages'][-1].content}")
                print()


if __name__ == "__main__":
    asyncio.run(main())
