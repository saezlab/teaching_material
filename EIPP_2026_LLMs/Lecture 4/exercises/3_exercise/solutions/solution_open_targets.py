"""
Solution 3 — Open Targets MCP
===============================
Connects to the Open Targets MCP server and uses an agent
to query drug and disease data for a gene.

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

    async with streamablehttp_client(server_url) as (read, write, _):
        async with ClientSession(read, write) as session:
            await session.initialize()

            tools = await load_mcp_tools(session)

            print(f"Found {len(tools)} tools:")
            for t in tools:
                print(f"  - {t.name}")
            print()

            llm = init_chat_model(
                model=MODEL_NAME,
                model_provider=MODEL_PROVIDER,
                base_url=BASE_URL,
                api_key=API_KEY,
                temperature=0.2,
            )

            agent = create_agent(model=llm, tools=tools)

            # Ask about diseases associated with TP53
            question = "What diseases are associated with TP53?"
            print(f"Question: {question}\n")

            result = await agent.ainvoke({
                "messages": [("user", question)]
            })

            print(f"Answer:\n{result['messages'][-1].content}")


if __name__ == "__main__":
    asyncio.run(main())
