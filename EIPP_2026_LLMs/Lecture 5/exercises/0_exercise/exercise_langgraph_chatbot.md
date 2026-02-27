# Exercise 0 — LangGraph Chatbot

## Goal

Build a simple chatbot using LangGraph's StateGraph. You will define the graph state, create a node that calls an LLM, wire the graph together, and run a chat loop.

## Prerequisites

- You have completed Lectures 2-4 exercises
- You have `langgraph` installed (if not: `uv add langgraph`)

## Setup

### 1. Create the exercise folder and copy the starter script

> [!TIP] Terminal
> ```bash
> mkdir -p lecture_5/0_langgraph_chatbot
> cp "../Lecture 5/exercises/0_exercise/exercise_langgraph_chatbot.py" lecture_5/0_langgraph_chatbot/
> ```

> [!NOTE]
> You can also copy the file using your file explorer. Just drag the Python file from `Lecture 5/exercises/0_exercise/` into `playground/lecture_5/0_langgraph_chatbot/`.

All commands below assume you are in the `playground/` folder.

---

## The Exercise

The starter script already provides:
- All necessary imports (`StateGraph`, `add_messages`, `init_chat_model`, etc.)
- A working chat loop that invokes the graph and prints responses

**Your tasks** (edit `lecture_5/0_langgraph_chatbot/exercise_langgraph_chatbot.py`):

1. **Define the `State` class** as a `TypedDict` with a single key `messages` using `Annotated[list, add_messages]`. This tells LangGraph to append new messages instead of overwriting the list.

2. **Implement the `call_llm` function** that:
   - Creates an LLM using `init_chat_model(model=MODEL_NAME, model_provider=MODEL_PROVIDER, base_url=BASE_URL, api_key=API_KEY)`
   - Invokes it with `state["messages"]`
   - Returns `{"messages": [ai_msg]}`

3. **Build the graph**:
   - Create a `StateGraph(State)`
   - Add a node called `"llm"` using your `call_llm` function
   - Set `"llm"` as the entry point
   - Add an edge from `"llm"` to `END`
   - Compile the graph and store it in the variable `graph`

> [!TIP] Terminal
> ```bash
> uv run lecture_5/0_langgraph_chatbot/exercise_langgraph_chatbot.py
> ```

**Expected outcome:** A working chatbot that remembers your conversation. Type messages and get responses. Type `quit` to exit.

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
    └── 0_langgraph_chatbot/
        └── exercise_langgraph_chatbot.py
```
