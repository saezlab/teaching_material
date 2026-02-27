# Exercise 1 — ReAct Agent with Biology Tools

## Goal

Build a ReAct agent using LangGraph that uses biology-themed tools to answer questions about DNA sequences. The agent will decide when to call tools, execute them, and use the results to produce a final answer.

## Prerequisites

- You have completed Lectures 2-4 exercises and Exercise 0 from this lecture
- You have `langgraph` installed (if not: `uv add langgraph`)

## Setup

### 1. Create the exercise folder and copy the starter script

> [!TIP] Terminal
> ```bash
> mkdir -p lecture_5/1_react_agent
> cp "../Lecture 5/exercises/1_exercise/exercise_react_agent.py" lecture_5/1_react_agent/
> ```

> [!NOTE]
> You can also copy the file using your file explorer. Just drag the Python file from `Lecture 5/exercises/1_exercise/` into `playground/lecture_5/1_react_agent/`.

All commands below assume you are in the `playground/` folder.

---

## The Exercise

The starter script already provides:
- Two biology tools: `gc_content` and `reverse_complement` (do not edit these)
- A `TOOLS` list containing both tools
- A main block that runs three questions through the agent

**Your tasks** (edit `lecture_5/1_react_agent/exercise_react_agent.py`):

1. **Define the `AgentState`** as a `TypedDict` with a single key `messages` using `Annotated[list, add_messages]`

2. **Implement the `agent_node` function** that:
   - Creates an LLM using `init_chat_model(model=MODEL_NAME, model_provider=MODEL_PROVIDER, base_url=BASE_URL, api_key=API_KEY)`
   - Binds the tools to the LLM: `llm_with_tools = llm.bind_tools(TOOLS)`
   - Invokes the LLM with `state["messages"]`
   - Returns `{"messages": [ai_msg]}`

3. **Implement the `route_after_agent` function** that:
   - Gets the last message from `state["messages"]`
   - Checks if it has tool calls using `getattr(last, "tool_calls", [])`
   - Returns `"tools"` if there are tool calls, `"end"` otherwise

4. **Build the ReAct graph**:
   - Create a `StateGraph(AgentState)`
   - Add node `"agent"` with your `agent_node` function
   - Add node `"tools"` with `ToolNode(TOOLS)`
   - Set `"agent"` as the entry point
   - Add conditional edges from `"agent"` using `route_after_agent`:
     ```python
     workflow.add_conditional_edges(
         "agent", route_after_agent, {"tools": "tools", "end": END}
     )
     ```
   - Add an edge from `"tools"` back to `"agent"` (this creates the loop!)
   - Compile the graph and store it in the variable `graph`

> [!TIP] Terminal
> ```bash
> uv run lecture_5/1_react_agent/exercise_react_agent.py
> ```

**Expected outcome:** The agent uses the biology tools to compute GC content and reverse complements. The third question requires the agent to call both tools in sequence.

---

## Bonus Challenge

Try adding your own question to the `questions` list that requires **both** tools, for example:
> "Is the GC content of ATGCGATCG higher than the GC content of its reverse complement?"

(Hint: the answer should be the same, since a sequence and its reverse complement always have equal GC content!)

---

## What your playground should look like now

```
playground/
├── .env
├── .venv/
├── pyproject.toml
├── lecture_2/
│   └── ...
├── lecture_3/
│   └── ...
├── lecture_4/
│   └── ...
└── lecture_5/
    ├── 0_langgraph_chatbot/
    │   └── exercise_langgraph_chatbot.py
    └── 1_react_agent/
        └── exercise_react_agent.py
```
