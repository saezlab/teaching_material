# Lecture 5: Agents and LangGraph

## Overview

This lecture introduces how to build agentic systems with Large Language Models (LLMs) using the ReAct pattern and LangGraph. You will start with core concepts (what an agent is, why ReAct helps), then learn the LangGraph primitives (state, nodes, edges, conditional routing, tool nodes), and finally implement a minimal ReAct loop end‑to‑end.

## Dependencies

> [!IMPORTANT]
> This lecture requires `langgraph`. If you completed Lecture 4 Exercise 2 (MCP), you already have it. Otherwise, run from your `playground/` folder:
> ```bash
> uv add langgraph
> ```

## Learning Objectives

By the end of this lecture, students will be able to:

- Understand what agents are and when to use them
- Explain the ReAct loop (Thought → Action → Observation) and its benefits
- Define typed graph state and implement nodes in LangGraph
- Connect nodes with edges and add conditional routing
- Add and execute tool nodes and invoke compiled graphs
- Replicate ReAct in LangGraph on a simple task

## Lecture Materials Index

### 1. Introduction to Agents and ReAct
**File:** [`0_agents_and_react.md`](0_agents_and_react.md)

Foundations and mental model: what agents are, why we use tools, and how ReAct structures reasoning and acting. Includes a diagram and the original paper citation.

**Key Topics:**
- What is an agent (LLM as a decision maker with tools)
- ReAct traces: Thought → Action → Observation → … → Final Answer
- When to use ReAct; common pitfalls (over‑tooling, infinite loops)
- Paper reference and conceptual grounding

---

### 2. LangGraph Basics
**File:** [`1_langgraph_basics.md`](1_langgraph_basics.md)

Hands‑on with LangGraph: define state, nodes, and edges; add conditional routes and a tool node; compile and invoke the graph.

**Key Topics:**
- `StateGraph`, typed state, message lists
- Nodes as Python callables returning partial state
- Edges and `add_conditional_edges(...)`
- Tool nodes (`ToolNode`) and execution of tool calls
- Compiling and invoking graphs

---

### 3. ReAct in LangGraph (End‑to‑End)
**File:** [`2_react_in_langgraph.md`](2_react_in_langgraph.md)

Put it together: bind tools to an LLM, route to a tool node when tool calls are present, loop back to the agent for the next Thought, and stop when no more tools are needed.

**Key Topics:**
- Binding tools to the model; inspecting `tool_calls`
- Routing after the agent node; returning to the agent after tools
- Stopping criteria; extracting the final answer
- Tracing and debugging the loop

---

## Prerequisites

- Completion of Lecture 2 (LLM APIs and LangChain), Lecture 3 (Structured Outputs) and Lecture 4 (Tool Use and MCP)
- Python 3.11+
- Basic familiarity with JSON, type hints, and virtual environments
- LangChain and LangGraph installed; provider API key configured to run LLM calls

## Getting Started

1. **Concepts First:** Read [`0_agents_and_react.md`](0_agents_and_react.md) and study the ReAct diagram and citation.
2. **LangGraph Basics:** Work through [`1_langgraph_basics.md`](1_langgraph_basics.md) to define state, nodes, edges, conditional routing, and a tool node.
3. **ReAct Loop:** Finish with [`2_react_in_langgraph.md`](2_react_in_langgraph.md) to run a minimal ReAct implementation and inspect its traces.
4. **Exercises:** Complete the exercises below to practice what you learned.

## Exercises

### Exercise 0 — LangGraph Chatbot
**File:** [`exercises/0_exercise/exercise_langgraph_chatbot.md`](exercises/0_exercise/exercise_langgraph_chatbot.md)

Build a simple chatbot using LangGraph's StateGraph. You will define the graph state, create a node that calls an LLM, wire the graph together, and run a chat loop.

### Exercise 1 — ReAct Agent with Biology Tools
**File:** [`exercises/1_exercise/exercise_react_agent.md`](exercises/1_exercise/exercise_react_agent.md)

Build a ReAct agent that uses biology-themed tools (GC content, reverse complement) to answer questions about DNA sequences. The agent decides when to call tools, executes them, and uses the results to produce a final answer.

## Practical Applications

- Domain Workflows: Orchestrate multi‑step tasks (search → read → compute → summarize).
- Fresh Knowledge Access: Call APIs/tools to retrieve up‑to‑date information.
- Operational Reliability: Keep control flow explicit and traces auditable.
- Portability: Reuse tools and graph patterns across providers and tasks.

## Notes on Agent/LangGraph Capabilities

Different providers support tool calls and message formats differently. LangGraph standardizes control flow while letting you bind provider‑specific models. Prefer typed state, explicit routing, and narrow tools for reliability and debuggability.

## Course Context

This lecture is part of the **EIPP 2026 LLMs** bootcamp. It builds on structured outputs (Lecture 3) and tool use (Lecture 4), and prepares you to implement robust agentic workflows.

## Additional Resources

- [ReAct Paper (arXiv)](https://arxiv.org/abs/2210.03629)
- [LangGraph Documentation](https://langchain-ai.github.io/langgraph/)
- [LangChain Tools Documentation](https://python.langchain.com/docs/how_to/#tools)
