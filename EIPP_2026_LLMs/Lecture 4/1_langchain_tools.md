# Using LangChain Tools

> [!NOTE]
> No new dependencies are needed for this section. The tools we define here use only pure Python and the LangChain packages you already have installed (`langchain`, `langchain-google-genai`, `langchain-openai`, `python-dotenv`).

## Table of Contents

1. [Using LangChain Tools](#using-langchain-tools)
2. [Define a tool](#define-a-tool)
3. [Bind tools to an LLM](#bind-tools-to-an-llm)
4. [Call the model and execute the tool call](#call-the-model-and-execute-the-tool-call)
5. [What we achieved](#what-we-achieved)

## Using LangChain Tools

In the previous section, we built a manual tool-calling system from scratch: we wrote helper functions to describe tools, crafted a system prompt with JSON schemas, and parsed the model's raw text output ourselves. That approach works, but it is tedious and error-prone.

LangChain provides a much cleaner way to do this. With a single `@tool` decorator, LangChain automatically generates the tool's name, description, and argument schema from your function. With `.bind_tools()`, the model knows exactly which tools are available and how to call them. And the model returns structured `tool_calls` that you can execute directly -- no manual JSON parsing needed.

In this guide, we'll walk through three key steps:

- **Define a tool** -- turn a plain Python function into a LangChain tool
- **Bind tools to an LLM** -- tell the model which tools it can use
- **Call the model and execute the tool call** -- let the model decide when to call a tool, run it, and return the result

We'll use two simple biology tools: one that calculates the **GC content** of a DNA sequence, and one that returns the **sequence length**. These are pure Python functions that any biologist will recognize immediately.

---

## Define a tool

In LangChain, a Python function becomes a tool by decorating it with `@tool`. The function's **docstring** tells the model what the tool does, and the **type hints** tell it what arguments to provide. This is all the model needs to decide when and how to call the tool.

Let's define two simple tools for working with DNA sequences:

```python
from langchain_core.tools import tool

@tool
def gc_content(sequence: str) -> float:
    """Calculate the GC content of a DNA sequence.

    GC content is the proportion of bases in a DNA sequence that are
    either guanine (G) or cytosine (C). It is a fundamental property
    used in molecular biology to characterize sequences.

    Args:
        sequence: A DNA sequence string containing A, T, G, and C bases.

    Returns:
        The GC content as a proportion between 0.0 and 1.0.
    """
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence)


@tool
def sequence_length(sequence: str) -> int:
    """Return the length of a DNA sequence.

    Counts the total number of bases (nucleotides) in the sequence.

    Args:
        sequence: A DNA sequence string containing A, T, G, and C bases.

    Returns:
        The number of bases in the sequence.
    """
    return len(sequence)
```

That's it. The `@tool` decorator turns each function into a LangChain tool object. LangChain reads the function name, docstring, and type hints to build a complete description that the model can understand.

You can inspect what LangChain generated:

```python
# The tool's name (taken from the function name)
print(gc_content.name)
# Output: gc_content

# The full argument schema (automatically generated from type hints)
print(gc_content.args_schema.model_json_schema())
# Output: {'description': 'Calculate the GC content of a DNA sequence. ...',
#          'properties': {'sequence': {'description': 'A DNA sequence string ...', 'type': 'string'}},
#          'required': ['sequence'], 'title': 'gc_contentSchema', 'type': 'object'}
```

Notice how the docstring and parameter descriptions end up in the schema. This is what the model reads to understand the tool. **Good docstrings lead to better tool use.**

You can also call the tools directly to verify they work:

```python
print(gc_content.invoke({"sequence": "ATGCGATCG"}))
# Output: 0.5555555555555556

print(sequence_length.invoke({"sequence": "ATGCGATCG"}))
# Output: 9
```

---

## Bind tools to an LLM

Now we need to tell the model about our tools. In LangChain, this is done with `.bind_tools()`. This attaches the tool schemas to the model so it knows what tools are available and how to call them.

```python
from dotenv import load_dotenv
load_dotenv()

from langchain.chat_models import init_chat_model

# --- Configuration ---
# MODEL_NAME = "gemini-3-flash-preview"
# MODEL_PROVIDER = "google_genai"
MODEL_NAME = "my-model"
MODEL_PROVIDER = "openai"
BASE_URL = "https://myung-courageous-gino.ngrok-free.dev/v1"
API_KEY = "not-needed"
llm = init_chat_model(model=MODEL_NAME, model_provider=MODEL_PROVIDER, base_url=BASE_URL, api_key=API_KEY, temperature=0.2)

# Bind both tools to the model
llm_with_tools = llm.bind_tools([gc_content, sequence_length])
```

That's the entire setup. The `llm_with_tools` object behaves just like the regular `llm`, but now the model knows it can call `gc_content` and `sequence_length` when it needs to.

A few things to note:

- You can bind as many tools as you want by adding them to the list.
- The model does **not** have to call a tool. If you ask it a question that doesn't need a tool (like "What is GC content?"), it will just answer normally.
- The model **decides on its own** when to use a tool based on the question and the tool descriptions.

---

## Call the model and execute the tool call

Now comes the key part. When we ask the model a question that requires computation, it will **not** try to compute the answer itself. Instead, it will emit a **tool call** -- a structured request saying "please run this tool with these arguments." We then execute the tool and (optionally) send the result back to the model so it can produce a final answer.

### Step 1: Ask the model a question

```python
question = "What is the GC content of the sequence ATGCGATCGATCGATCGATCG?"

ai_msg = llm_with_tools.invoke(question)
```

### Step 2: Inspect the tool calls

The model's response is an `AIMessage`. If the model decided to use a tool, the message will contain one or more `tool_calls`:

```python
print(ai_msg.tool_calls)
# Output: [{'name': 'gc_content',
#           'args': {'sequence': 'ATGCGATCGATCGATCGATCG'},
#           'id': '...',
#           'type': 'tool_call'}]
```

Each tool call has:
- `name` -- which tool the model wants to call
- `args` -- the arguments it wants to pass
- `id` -- a unique ID that links the tool call to its result

If the model thinks it can answer without a tool, `tool_calls` will be an empty list. You should check for this before trying to execute tool calls.

### Step 3: Execute the tool

We look up the tool by name and run it with the arguments the model provided:

```python
# Build a lookup dictionary so we can find tools by name
tools_by_name = {
    gc_content.name: gc_content,
    sequence_length.name: sequence_length,
}

# Execute each tool call and store the results
results = []
for tool_call in ai_msg.tool_calls:
    tool = tools_by_name[tool_call["name"]]
    result = tool.invoke(tool_call["args"])
    results.append((tool_call, result))
    print(f"Tool: {tool_call['name']}, Result: {result}")
```

At this point, we already have the answer (the GC content value). But often we want the model to **interpret** the result and produce a human-friendly response. To do this, we send the tool result back to the model.

### Step 4: Return the result to the model for a final answer

We rebuild the full conversation -- the original question, the model's tool call, and the tool's result -- and send it back to the model. The model then uses the tool output to write a final answer.

A `ToolMessage` tells the model "here is the result from the tool you asked to call."
The `tool_call_id` links each result to the specific tool call that requested it --
this is required so the model knows which result belongs to which request.

```python
from langchain_core.messages import HumanMessage, ToolMessage

# Build messages with stored results
messages = [HumanMessage(content=question)]  # The original question
messages.append(ai_msg)                       # The model's tool call

# Add the tool result for each tool call
for tool_call, result in results:
    messages.append(
        ToolMessage(
            content=str(result),
            tool_call_id=tool_call["id"],
        )
    )

# Use the base llm (not llm_with_tools) so the model writes a text answer
# instead of potentially trying to call more tools
final_answer = llm.invoke(messages)
print(final_answer.content)
```

The model will now produce something like:

> The GC content of the sequence ATGCGATCGATCGATCGATCG is approximately 0.524, meaning about 52.4% of the bases are guanine (G) or cytosine (C).

This is the complete tool-calling loop:

1. **You ask** a question
2. **The model decides** to call a tool and emits a structured tool call
3. **You execute** the tool with the provided arguments
4. **You send the result back** to the model
5. **The model produces** a final human-readable answer

---

## What we achieved

1. We defined two simple biology tools (`gc_content` and `sequence_length`) using the `@tool` decorator, with clear docstrings that help the model understand when and how to use them.
2. We bound the tools to an LLM with `.bind_tools()`, giving the model access to our tools without any manual prompt engineering or JSON schema crafting.
3. We invoked the model, captured its structured tool call, executed the tool, and returned the result so the model could produce a final narrative answer.

Compare this to the manual approach from the previous section: no more hand-written tool descriptions, no more JSON parsing, no more hoping the model follows your output format. LangChain handles all of that automatically. This pattern -- **define tool, bind to LLM, execute tool calls** -- is the foundation for everything we'll build in the rest of this lecture.

## Exercise

Time to put this into practice! See [exercise_langchain_tools.md](exercises/1_exercise/exercise_langchain_tools.md) — you will define biology-themed tools, bind them to an LLM, and execute tool calls.
