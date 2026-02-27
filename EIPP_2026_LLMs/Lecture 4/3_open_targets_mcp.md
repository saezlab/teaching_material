# Querying Open Targets with MCP

## Table of Contents

1. [What is Open Targets?](#what-is-open-targets)
2. [Connect to the Open Targets MCP server](#connect-to-the-open-targets-mcp-server)
3. [Build an agent that queries Open Targets](#build-an-agent-that-queries-open-targets)
4. [What happens under the hood](#what-happens-under-the-hood)
5. [Verify on the platform](#verify-on-the-platform)

> [!NOTE]
> This section requires the MCP packages installed in the previous section (`uv add mcp langchain-mcp-adapters langgraph`). If you have already done this, no action needed.

## What is Open Targets?

The [Open Targets Platform](https://platform.opentargets.org) is a public resource for systematic drug target identification and prioritisation. It integrates genetics, genomics, transcriptomics, and drug data into a single portal where you can look up:

- **Target–disease associations** (e.g. "which diseases is BRAF linked to?")
- **Known drugs** for a target (e.g. "what drugs target BRAF?")
- **Genetic evidence**, pathways, and more

Normally you would search Open Targets by visiting their website. The MCP server lets your Python code do the same thing automatically -- and an agent can chain multiple queries to answer complex questions.

Open Targets exposes a remote **MCP server** at `https://mcp.platform.opentargets.org/mcp`. This means we can connect to it with the same `streamablehttp_client` pattern we used for BioContext.AI — and let an agent query the platform programmatically.

## Connect to the Open Targets MCP server

First, let's connect and inspect which tools the server exposes.

```python
import asyncio
from mcp import ClientSession
from mcp.client.streamable_http import streamablehttp_client
from langchain_mcp_adapters.tools import load_mcp_tools

async def main():
    server_url = "https://mcp.platform.opentargets.org/mcp"

    async with streamablehttp_client(server_url) as (read, write, _):
        async with ClientSession(read, write) as session:
            await session.initialize()

            tools = await load_mcp_tools(session)

            print(f"Found {len(tools)} tools:")
            for t in tools:
                print(f"  - {t.name}: {t.description[:80]}...")

if __name__ == "__main__":
    asyncio.run(main())
```

You should see several tools related to searching targets, diseases, drugs, and associations. The exact set may change as the platform evolves, but typical tools include search, target lookup, drug queries, and association retrieval.

> [!NOTE]
> If you see a connection error, check your internet connection and try again
> in a few minutes.

## Build an agent that queries Open Targets

Now let's combine the MCP tools with `create_agent` to ask a real biology question: **"What drugs target BRAF?"**

```python
import asyncio
from mcp import ClientSession
from mcp.client.streamable_http import streamablehttp_client
from langchain_mcp_adapters.tools import load_mcp_tools
from langchain.agents import create_agent
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

async def main():
    server_url = "https://mcp.platform.opentargets.org/mcp"

    async with streamablehttp_client(server_url) as (read, write, _):
        async with ClientSession(read, write) as session:
            await session.initialize()

            tools = await load_mcp_tools(session)

            llm = init_chat_model(
                model=MODEL_NAME,
                model_provider=MODEL_PROVIDER,
                base_url=BASE_URL,
                api_key=API_KEY,
                temperature=0.2,
            )

            agent = create_agent(model=llm, tools=tools)

            result = await agent.ainvoke({
                "messages": [("user", "What drugs target BRAF?")]
            })

            print(result["messages"][-1].content)

if __name__ == "__main__":
    asyncio.run(main())
```

## What happens under the hood

When the agent processes "What drugs target BRAF?", it typically follows a multi-step reasoning flow:

1. **Search for BRAF** — The agent calls a search tool to look up "BRAF" and resolves it to a specific identifier (Ensembl gene ID `ENSG00000157764`).
2. **Query drug data** — Using the resolved ID, the agent calls a tool that retrieves known drugs targeting this gene (e.g. vemurafenib, dabrafenib, encorafenib, sorafenib).
3. **Summarise** — The agent formats the results into a human-readable answer.

The agent decides which tools to call and in what order — you just ask the question. This is the power of combining MCP tools with an agent loop: the model can chain multiple tool calls to answer complex queries.

> **Note:** The exact tool names and call sequence depend on what the Open Targets MCP server exposes. The agent adapts its strategy based on the available tools.

## Verify on the platform

To check that the agent's answer is correct, visit the Open Targets Platform directly:

**https://platform.opentargets.org/target/ENSG00000157764**

Navigate to the **Known drugs** section. You should see the same drugs the agent reported (vemurafenib, dabrafenib, etc.), confirming that the MCP tools are returning real, up-to-date data from the platform.

This is a good habit when working with any external data source: always cross-check a sample of the results manually.

## Exercise

Time to try it yourself! See [exercise_open_targets.md](exercises/3_exercise/exercise_open_targets.md) — you will connect to the Open Targets MCP server and query it with your own biology question.
