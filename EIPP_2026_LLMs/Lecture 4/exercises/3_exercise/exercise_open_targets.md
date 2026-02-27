# Exercise 3 — Querying Open Targets with MCP

## Goal

Connect to the remote Open Targets MCP server and use an agent to answer a biology question about a gene of your choice. Then cross-check the answer on the Open Targets Platform website.

## Prerequisites

- You have read [`3_open_targets_mcp.md`](../../3_open_targets_mcp.md)
- You have completed exercises 0–2 of this lecture
- MCP packages are already installed (`mcp`, `langchain-mcp-adapters`, `langgraph`)

## Setup

### 1. Create the exercise folder and copy the starter file

> [!TIP] Terminal
> ```bash
> mkdir -p lecture_4/3_open_targets
> cp "../Lecture 4/exercises/3_exercise/exercise_open_targets.py" lecture_4/3_open_targets/
> ```

All commands below assume you are in the `playground/` folder.

---

## The Exercise

You are given a starter script that connects to the Open Targets MCP server. Your job is to complete it and ask a biology question.

### Part A: Connect and inspect tools

**Your tasks** (edit `lecture_4/3_open_targets/exercise_open_targets.py`):

1. **Complete the MCP connection** using `streamablehttp_client` to connect to `https://mcp.platform.opentargets.org/mcp`
2. **Load the tools** with `load_mcp_tools(session)` and print their names to see what's available

### Part B: Build the agent and ask a question

3. **Create an agent** with `create_agent(model=llm, tools=tools)`
4. **Ask a question about a gene** — choose any gene that interests you. Some ideas:
   - "What diseases are associated with TP53?"
   - "What drugs target EGFR?"
   - "What is known about the BRCA1 gene?"
5. **Run the script** and read the agent's answer

> [!NOTE]
> If the Open Targets MCP server is unreachable, wait a few minutes and try again.
> You can check your internet connection with `curl https://mcp.platform.opentargets.org/mcp`.

> [!TIP] Terminal
> ```bash
> uv run lecture_4/3_open_targets/exercise_open_targets.py
> ```

### Part C: Cross-check on the website

6. **Go to the Open Targets Platform** at https://platform.opentargets.org
7. **Search for the same gene** you asked about
8. **Compare** the agent's answer with the website. Look for: the same drug names, the same disease associations, or the same gene descriptions. Minor differences in wording are expected, but the key facts should match.

---

## Bonus Challenge

Try a question that requires the agent to chain multiple tool calls. For example:
> "Which drugs that target BRAF are approved for melanoma?"

Watch the agent's reasoning — does it search for BRAF first, then filter by disease?

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
│   ├── 2_mcp/
│   │   ├── math_server.py
│   │   └── exercise_mcp_client.py
│   └── 3_open_targets/
│       └── exercise_open_targets.py
```
