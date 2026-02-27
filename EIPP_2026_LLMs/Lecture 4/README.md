# Lecture 4: Tool Use and Model Context Protocol (MCP)

## Overview

This lecture introduces how to supercharge Large Language Models (LLMs) with external tools and the Model Context Protocol (MCP). You will start with first principles (having a model emit structured calls to simple Python functions), then move to LangChain’s native tool interface, and finally learn how MCP standardizes tool definitions so they can be reused across models and agent frameworks.

## Learning Objectives

By the end of this lecture, students will be able to:

- Understand why tools are needed to overcome LLM limitations (knowledge cutoff, precision)
- Prompt a model to emit strict JSON for tool invocation and safely parse it
- Define Python tools with LangChain’s `@tool` decorator and bind them to an LLM
- Execute tool calls emitted by an LLM and return results for final synthesis
- Explain MCP concepts (typed tools, transports, client/server separation)
- Build a minimal MCP server (FastMCP) and call it from a LangChain agent (LangGraph powers the agent loop internally and will be covered in depth in Lecture 5)
- Connect to a remote MCP server, filter the tools you need, and execute them
- Query real biomedical databases (Open Targets) through remote MCP servers

## Lecture Materials Index

### 1. Introduction to Tools
**📁 File:** [`0_introduction_to_tools.md`](0_introduction_to_tools.md)

Ground-up approach to tool use: describe Python functions, expose their signatures to the model, and require strict JSON outputs that match the tool’s parameters.

**Key Topics:**
- Why tools: access to fresh knowledge and precise computation
- Describing functions (docstrings, signatures) to the model
- Building a system prompt with a tool registry and JSON schema
- Parsing model outputs into executable calls and returning results

---

### 2. Using LangChain Tools
**📁 File:** [`1_langchain_tools.md`](1_langchain_tools.md)

Define and run tools the idiomatic LangChain way, using a biology-themed example (GC content analysis) to demonstrate definition, binding, and execution.

**Key Topics:**
- `@tool` decorator: schema from signature + docstring
- `bind_tools(...)` to declare available tools to the LLM
- Handling `tool_calls` and executing them safely
- Returning tool outputs and prompting the model for final summaries

---

### 3. Introduction to Model Context Protocol (MCP)
**📁 File:** [`2_intro_to_mcp.md`](2_intro_to_mcp.md)

Make tools portable and vendor‑neutral. Create a small FastMCP server that exposes `add` and `multiply`, then call it from a LangChain + LangGraph agent. Also learn how to query a remote MCP server and invoke domain tools.

**Key Topics:**
- What MCP is and why it matters (typed I/O, portability, composability)
- FastMCP server exposing typed tools over `stdio`
- Loading MCP tools into LangChain and using an agent to call them
- Connecting to a remote MCP server and filtering tools by name

---

### 4. Querying Open Targets with MCP
**📁 File:** [`3_open_targets_mcp.md`](3_open_targets_mcp.md)

Connect to a real biomedical database through a remote MCP server. Build an agent that queries the Open Targets Platform to answer questions about drug targets and disease associations.

**Key Topics:**
- The Open Targets Platform as a drug-target knowledge base
- Connecting to a remote MCP server with `streamablehttp_client`
- Building an agent with `create_agent` to query biological data
- Cross-checking agent answers against the platform website

---

## Prerequisites

- Completion of Lecture 2 (LLM APIs and LangChain) and Lecture 3 (Structured Outputs)
- Python 3.11 with UV (set up in Lecture 2)
- `.env` file in `playground/` (set up in Lecture 2)
- Basic familiarity with JSON, type hints, and virtual environments
- MCP packages will be installed during the lecture (see Getting Started)

## Getting Started

1. **Install MCP dependencies:** From the `playground/` directory, run `uv add mcp langchain-mcp-adapters langgraph` to install the packages needed for the MCP section.
2. **Introduction to Tools:** Read [`0_introduction_to_tools.md`](0_introduction_to_tools.md) — then do **Exercise 0**.
3. **LangChain Tools:** Read [`1_langchain_tools.md`](1_langchain_tools.md) — then do **Exercise 1**.
4. **MCP Basics:** Read [`2_intro_to_mcp.md`](2_intro_to_mcp.md) — then do **Exercise 2**.
5. **Open Targets MCP:** Read [`3_open_targets_mcp.md`](3_open_targets_mcp.md) — then do **Exercise 3**.

Each guide ends with a link to its corresponding exercise. Do the exercise right after the guide so the concepts are fresh.

> **Note:** The MCP exercises (2 and 3) use Python's `async`/`await` syntax. This is introduced in the MCP guide -- follow the patterns even if the syntax is unfamiliar.

## Exercises

- **Exercise 0 — Manual Tool Calling** ([`exercise_manual_tools.md`](exercises/0_exercise/exercise_manual_tools.md)): Extend the manual calculator with a divide tool and error handling.
- **Exercise 1 — LangChain Tools** ([`exercise_langchain_tools.md`](exercises/1_exercise/exercise_langchain_tools.md)): Define biology-themed tools, bind them to an LLM, and execute tool calls.
- **Exercise 2 — MCP Server and Client** ([`exercise_mcp.md`](exercises/2_exercise/exercise_mcp.md)): Extend a math MCP server with a new tool and build a client.
- **Exercise 3 — Open Targets MCP** ([`exercise_open_targets.md`](exercises/3_exercise/exercise_open_targets.md)): Query the Open Targets Platform with an agent and cross-check results.

## Practical Applications

- **Domain Workflows:** Let models call bioinformatics/R&D tools (enrichment, querying knowledge bases) with typed inputs/outputs.
- **Fresh Knowledge Access:** Bridge the knowledge cutoff by calling web/data tools.
- **Operational Reliability:** Keep business logic in Python; let the model orchestrate via tool calls.
- **Portability:** Define once with MCP, reuse across providers and agents.

## Notes on Tool/MCP Capabilities

Different providers offer varying support for tool calls and structured outputs. LangChain abstracts many differences, while MCP standardizes tool description and invocation across runtimes. Prefer typed schemas and explicit argument validation for reliability.

## Course Context

This lecture is part of the **EIPP 2026 LLMs bootcamp**. It builds on structured outputs (Lecture 3) and prepares you to integrate external computation and data sources safely in later lectures.

## Additional Resources

- [LangChain Tools Documentation](https://python.langchain.com/docs/how_to/#tools)
- [LangGraph ReAct Agents](https://langchain-ai.github.io/langgraph/)
- [Model Context Protocol (MCP) Docs](https://modelcontextprotocol.io)
- [MCP GitHub Organization](https://github.com/modelcontextprotocol)


