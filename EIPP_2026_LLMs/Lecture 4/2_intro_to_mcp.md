# Introduction to Model Context Protocol (MCP)

## Table of Contents

1. [Install the required packages](#install-the-required-packages)
2. [Model Context Protocol (MCP): a gentle introduction](#model-context-protocol-mcp-a-gentle-introduction)
3. [Why MCP?](#why-mcp)
4. [Create a minimal MCP server (FastMCP)](#create-a-minimal-mcp-server-fastmcp)
5. [Call the MCP tools from a LangChain + LangGraph agent](#call-the-mcp-tools-from-a-langchain--langgraph-agent)
6. [Call a remote MCP server](#call-a-remote-mcp-server)


## Install the required packages

Before writing any code, install the three packages we need:

> [!TIP] Terminal
> ```bash
> uv add mcp langchain-mcp-adapters langgraph
> ```

What each package does:
- **`mcp`** -- the Model Context Protocol SDK (server + client).
- **`langchain-mcp-adapters`** -- bridges MCP tools into LangChain so you can use them like any other LangChain tool.
- **`langgraph`** -- the runtime that powers the agent loop behind the scenes.

## Model Context Protocol (MCP): a gentle introduction

So far, tools live inside your Python script. But what if you want to share tools between projects, or use tools written by someone else? That's what MCP solves.

MCP is a simple, open protocol that lets AI models securely access tools and data sources through a standard interface. Instead of every app inventing its own way to “wire up” tools, MCP defines how tools are described, discovered, and called, as well as how results and errors are returned. The goal is portability (use the same tool with many model runtimes/agents), safety (typed inputs/outputs), and simplicity (one mental model, multiple transports such as stdio/HTTP/SSE).

- **Docs**: [modelcontextprotocol.io](https://modelcontextprotocol.io)
- **Spec and examples**: [github.com/modelcontextprotocol](https://github.com/modelcontextprotocol)

### Why MCP?

- **Vendor‑neutral tools**: you can reuse the same tool with different LLM providers or agent frameworks.
- **Typed I/O**: tools come with schemas so models know how to call them.
- **Composability**: mix MCP tools with other LangChain tools in one agent.

A **transport** is how the client and server communicate -- for example, over standard input/output (stdio) for local processes, or over HTTP for remote servers.

In this short guide, we’ll build a tiny math tool server with MCP and then call it from a LangChain + LangGraph agent. The narrative is intentionally simple: we want the model to answer “What is 3 + 5?” by deciding to call an MCP tool rather than guessing.

## Create a minimal MCP server (FastMCP)

We’ll expose two tools, `add` and `multiply`, over MCP using the stdio transport. Save this as `math_server.py` and run it with Python.

```python
# math_server.py
from mcp.server.fastmcp import FastMCP

mcp = FastMCP("Math")

@mcp.tool()
def add(a: int, b: int) -> int:
    """Add two numbers"""
    return a + b

@mcp.tool()
def multiply(a: int, b: int) -> int:
    """Multiply two numbers"""
    return a * b

if __name__ == "__main__":
    # Expose the server over stdio (simple and reliable for local development)
    mcp.run(transport="stdio")
```

What this gives us:
- The server announces a capability set and a tool registry to clients.
- Each tool has a name, a short description, and typed parameters.
- Any MCP‑aware client (not just LangChain) can discover and call these tools.

## Call the MCP tools from a LangChain + LangGraph agent

> [!NOTE]
> **About async/await**
>
> MCP uses Python's asynchronous programming features because it needs to communicate
> with external servers (sending requests and waiting for responses). You will see three
> new keywords:
> - `async def` -- marks a function as asynchronous (it can "wait" for things)
> - `await` -- pauses execution until an operation completes (like waiting for a server response)
> - `async with` -- like regular `with`, but for async resources (connections, sessions)
>
> At the end of the script, `asyncio.run(main())` starts the async event loop and runs
> your `main()` function. You do not need to deeply understand async programming --
> just follow the pattern. All MCP code in this course uses the same structure.
>
> If you are running this code in a Jupyter notebook (not used in this course by default),
> you will need to install and apply `nest_asyncio` since notebooks already run an event loop:
> ```python
> import nest_asyncio
> nest_asyncio.apply()
> ```

> [!NOTE]
> **What is an agent?**
>
> An **agent** is a program that uses an LLM to decide which tools to call and in what
> order. It works in a loop:
> 1. **Reason** -- the model looks at the question and available tools
> 2. **Act** -- the model calls a tool with specific arguments
> 3. **Observe** -- the model reads the tool's result
> 4. **Repeat** -- if more information is needed, go back to step 1
>
> This loop is called **ReAct** (Reason + Act). The `create_agent` function from
> LangChain sets up this loop automatically -- you just provide a model and tools.
> In the previous guides (and Exercises 0-1), you did this loop manually (parse tool calls, execute, return results).
> From now on, `create_agent` handles it for you. The logic is the same, just automated.

Below we connect to the server via stdio, load the exposed tools into LangChain using `langchain_mcp_adapters`, and create an agent with `create_agent` from `langchain.agents`. This is LangChain's high-level agent API -- it wires up a ReAct-style loop under the hood, so you just pass a model and tools. The agent will choose when to call `add`.

```python
# client_langchain_mcp.py
import asyncio
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from langchain_mcp_adapters.tools import load_mcp_tools
from langchain.agents import create_agent

# Uses the helper introduced earlier in these notes
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
    import os
    server_path = os.path.join(os.path.dirname(__file__), "math_server.py")
    server_params = StdioServerParameters(
        command="uv",
        args=["run", server_path],
    )

    async with stdio_client(server_params) as (read, write):
        async with ClientSession(read, write) as session:
            # Initialize the MCP connection and discover tools
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

            # Ask a question that should trigger a call to the `add` tool
            result = await agent.ainvoke({
                "messages": [("user", "What is 3 + 5?")]
            })

            # Show the final answer the agent produced
            print(result["messages"][-1].content)

if __name__ == "__main__":
    asyncio.run(main())
```

### How this works (conceptually)

1. We run an MCP server that exports typed tools over stdio.
2. The client opens a session, initializes the protocol, and loads tools with `load_mcp_tools(...)`.
3. The agent reasons about the user request and selects one of the MCP tools (e.g., `add`) with appropriate arguments.
4. The tool executes in the server, returns a typed result, and the agent uses it to produce the final answer.

### Tips and troubleshooting

- Use an absolute path for `math_server.py` in `StdioServerParameters` -- the code above does this with `os.path.dirname(__file__)`.
- Using `uv run` ensures the server runs inside the correct virtual environment with all dependencies available.
- With stdio transport, the client starts and stops the server automatically -- you don't need to run the server separately.
- You can mix MCP tools with any other LangChain tools in the same agent.

That’s all you need to get started: define a tool once with MCP, reuse it across models and frameworks, and keep your integration surface clean and portable.

## Call a remote MCP server

In the previous section we saw how we can define a local MCP server. More precisely, we created a script that gets started and exposes two tools over MCP. In this section we will see how we can call a remote MCP server (e.g. a server running on a remote machine). Crucially, this step leverages a resource called **[BioContext.AI](https://biocontext.ai)** that we recently released in collaboration with other researchers. The server is hosted at `https://biocontext-kb.fastmcp.app/mcp`.

> [!NOTE]
> In this section we use the manual `bind_tools` approach (like Exercise 1) to show
> how remote MCP tools work step by step. In the next section (Open Targets), we use
> `create_agent` which handles the tool-calling loop automatically. Both approaches
> work -- the manual approach gives you more control, while `create_agent` is simpler.

The following code will allow us to query programmatically the STRING database through the MCP server and retrieve the interactors of TP53 in human.

```python
import asyncio
from mcp import ClientSession
from mcp.client.streamable_http import streamablehttp_client
from langchain_mcp_adapters.tools import load_mcp_tools
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
    mcp_server = "https://biocontext-kb.fastmcp.app/mcp"

    # The _ is a third value we don't need (session configuration)
    async with streamablehttp_client(mcp_server) as (read, write, _):
        async with ClientSession(read, write) as session:
            # Initialize the connection
            await session.initialize()

            tools = await load_mcp_tools(session)

            # We subset only on the tools from STRING
            string_tools = [tool for tool in tools if "string" in tool.name]

            # Define the LLM
            llm = init_chat_model(model=MODEL_NAME, model_provider=MODEL_PROVIDER, base_url=BASE_URL, api_key=API_KEY)

            # Bind the tools to the LLM
            llm_with_tools = llm.bind_tools(string_tools)

            # Invoke the LLM
            res = await llm_with_tools.ainvoke("What are the interactors of TP53 in human?")

            # Obtain the tool calls
            tool_calls = res.tool_calls

            # Execute the tool calls
            for call in tool_calls:
                name = call["name"]
                args = call.get("args", {})
                for tool in string_tools:
                    if tool.name == name:
                        result = await tool.ainvoke(args)
                        print(result)

if __name__ == "__main__":
    asyncio.run(main())
```

## Exercise

Time to put this into practice! See [exercise_mcp.md](exercises/2_exercise/exercise_mcp.md) — you will extend a math MCP server with a new tool and build a client that calls it.

