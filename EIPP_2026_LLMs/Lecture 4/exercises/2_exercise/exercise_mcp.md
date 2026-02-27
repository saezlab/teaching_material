# Exercise 2 — MCP Server and Client

## Goal

Practice building an MCP server with FastMCP, extending it with a new tool, and calling it from a LangChain + LangGraph agent.

## Prerequisites

- You have read [`2_intro_to_mcp.md`](../../2_intro_to_mcp.md)
- You have completed exercises 0 and 1 of this lecture

## Setup

### 1. Install the MCP packages (if not already done)

> [!TIP] Terminal
> ```bash
> uv add mcp langchain-mcp-adapters langgraph
> ```

### 2. Create the exercise folder and copy the starter files

> [!TIP] Terminal
> ```bash
> mkdir -p lecture_4/2_mcp
> cp "../Lecture 4/exercises/2_exercise/math_server.py" lecture_4/2_mcp/
> cp "../Lecture 4/exercises/2_exercise/exercise_mcp_client.py" lecture_4/2_mcp/
> ```

> [!NOTE]
> You can also copy the files using your file explorer. Just drag both Python files from `Lecture 4/exercises/2_exercise/` into `playground/lecture_4/2_mcp/`.

All commands below assume you are in the `playground/` folder.

---

## The Exercise

You are given a working MCP math server (`math_server.py`) with `add` and `multiply` tools. You need to extend it and build a client that uses it.

### Part A: Extend the server

**Your tasks** (edit `lecture_4/2_mcp/math_server.py`):

1. **Add a `subtract` tool** that subtracts two integers (returns `a - b`)
2. **Add a `power` tool** that raises `a` to the power of `b` (returns `a ** b`)
3. Both tools need the `@mcp.tool()` decorator and a clear docstring

### Part B: Build the client

**Your tasks** (edit `lecture_4/2_mcp/exercise_mcp_client.py`):

1. **Set up the server connection** using `StdioServerParameters` pointing to your `math_server.py` (`stdio` stands for "standard input/output" -- a way for two programs on the same machine to communicate by sending text back and forth)
2. **Load the MCP tools** into LangChain using `load_mcp_tools(session)`
3. **Create an agent** with `create_agent(model=llm, tools=tools)`
4. **Test with a question** that requires tool use, e.g., "What is 2 to the power of 10?"

> [!TIP] Terminal
> ```bash
> uv run lecture_4/2_mcp/exercise_mcp_client.py
> ```

**Expected outcome:** The agent discovers the available tools from the MCP server, selects the right tool (e.g., `power`), calls it with the correct arguments, and returns the result (1024). You may see `Processing request of type ...` messages in the output — these are normal debug output from the MCP server and can be ignored.

> [!TIP] Troubleshooting
> If the script hangs or crashes, check your `math_server.py` for syntax errors.
> You can test the server independently by running `uv run lecture_4/2_mcp/math_server.py`
> -- if it starts without errors, it is working correctly (press Ctrl+C to stop).
> The server won't print anything — that's normal. It's waiting for a client to connect via stdin.

---

## Bonus Challenge

Try asking a multi-step question like:
> "What is 3 + 5, and then raise that result to the power of 2?"

Does the agent chain the tool calls correctly?

---

## What your playground should look like now

```
playground/
├── .env
├── .venv/
├── pyproject.toml
├── lecture_4/
│   ├── 0_manual_tools/
│   │   └── exercise_manual_tools.py
│   ├── 1_langchain_tools/
│   │   └── exercise_langchain_tools.py
│   └── 2_mcp/
│       ├── math_server.py
│       └── exercise_mcp_client.py
```
