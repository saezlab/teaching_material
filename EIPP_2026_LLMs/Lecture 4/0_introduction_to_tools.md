# Introduction to Tools

> [!NOTE]
> This guide assumes you already have the playground set up from Lecture 2 with these dependencies installed: `langchain`, `langchain-google-genai`, `langchain-openai`, and `python-dotenv`. If not, go back to Lecture 2, Exercise 0 and get set up first.

## Table of Contents

1. [Introduction to Tools](#introduction-to-tools)
2. [A simple tool](#a-simple-tool)

## Introduction to Tools

Large Language Models (LLMs) become far more useful when combined with external tools. Tools allow them to overcome key limitations and extend their reasoning capabilities in practical ways:

* **Access to new knowledge**: LLMs have a knowledge cutoff, so tools let them fetch up-to-date or domain-specific information that isn’t in their training data.
* **Task offloading**: Some problems are not easily solved by text-based reasoning; dedicated tools (e.g., calculators, databases, APIs) can handle them more efficiently, much like humans rely on calculators for precise math.

## A simple tool

Every function that takes as input a set of arguments and returns a result can be seen as a tool. However, we must recall that LLMs interact with the world through text, so in order for an LLM to use a tool, we must be able to describe it in text. By using what we learned in the previous lecture, we can try to have a model use two tools: a tool to do additions and a tool to do multiplications. Let us first define a function that adds two numbers together.

```python
def add(a: float, b: float) -> float:
    '''Add two numbers together'''
    return a + b

def multiply(a: float, b: float) -> float:
    '''Multiply two numbers together'''
    return a * b
```

It is important to note how we defined the functions above. We added a docstring to the function that describes what the function does. This would be key to have the model understand what the function does and how to use it.

Before we go further, let's set up our imports and load our API key:

```python
from dotenv import load_dotenv
load_dotenv()

import json
import inspect
from typing import Callable, List, get_type_hints

from langchain.chat_models import init_chat_model
from langchain_core.prompts import ChatPromptTemplate
```

Now let us prepare a system prompt that will describe the tools to the model, making it aware of the tools it can use. In order to do so we need to make sure that the model writes a JSON object matching the signature of the function it wants to use.

The model can only read text, so we need to convert our Python function signatures into a textual description that the model can understand.

The four helper functions below do the following:
- `_param_spec(func)` -- reads a function's parameter names and types
- `build_tools_block(funcs)` -- creates a text description of all tools for the system prompt
- `build_json_schema(funcs)` -- creates a JSON schema the model must follow
- `build_system_prompt(funcs)` -- combines everything into a complete system prompt

The helper functions below handle the plumbing of building tool descriptions and JSON schemas. **You do not need to memorize this code** — read through it to understand the idea, and in the next section we'll see how LangChain does all of this for you automatically.

```python
def _param_spec(func: Callable) -> List[dict]:
    sig = inspect.signature(func)
    hints = get_type_hints(func)
    specs = []
    for name, p in sig.parameters.items():
        typ = hints.get(name, str)
        specs.append({
            "name": name,
            "type": getattr(typ, "__name__", str(typ)),
            "required": p.default is inspect._empty,
            "default": None if p.default is inspect._empty else p.default,
        })
    return specs

def build_tools_block(funcs: List[Callable]) -> str:
    lines = []
    for f in funcs:
        params = _param_spec(f)
        lines.append(f"- {f.__name__}: {f.__doc__ or ''}".strip())
        for p in params:
            default_str = "" if p["required"] else f" (default={p['default']})"
            lines.append(f"    * {p['name']}: {p['type']}{default_str}")
    return "\n".join(lines)

def build_json_schema(funcs: List[Callable]) -> str:
    names = [f.__name__ for f in funcs]
    return json.dumps({
        "type": "object",
        "properties": {
            "tool": {"type": ["string", "null"], "enum": names + [None]},
            "args": {"type": "object"}
        },
        "required": ["tool", "args"],
        "additionalProperties": False
    }, indent=2)

def build_system_prompt(funcs: List[Callable]) -> str:
    tools_block = build_tools_block(funcs)
    schema = build_json_schema(funcs)
    examples = [
        {"tool": "add", "args": {"a": 1, "b": 2}},
        {"tool": "multiply", "args": {"a": 3, "b": 4}},
    ]
    examples_str = "\n".join(json.dumps(e, indent=2) for e in examples)
    return f"""You can call Python tools by emitting exactly ONE JSON object with this shape:

- tool: one of the tool names below (or null if no tool is needed)
- args: an object whose keys match the tool's parameters exactly

TOOLS AVAILABLE:
{tools_block}

OUTPUT FORMAT (STRICT):
Return ONLY a JSON object, no prose, no code fences.

Call object schema (guidance):
{schema}

EXAMPLES (strict JSON, no extra text):
{examples_str}

If you do NOT need any tool, append to the end of the response:
{{"tool": null, "args": {{}}}}
"""
```

We can try now to see one example of system prompt.

```python
funcs = [add, multiply]
system_prompt = build_system_prompt(funcs)
print(system_prompt)
```

We can now invoke the model with the system prompt we just prepared.

```python
# --- Configuration ---
# MODEL_NAME = "gemini-3-flash-preview"
# MODEL_PROVIDER = "google_genai"
MODEL_NAME = "my-model"
MODEL_PROVIDER = "openai"
BASE_URL = "https://myung-courageous-gino.ngrok-free.dev/v1"
API_KEY = "not-needed"

prompt = ChatPromptTemplate.from_messages([
    ("system", "{system_prompt}"),
    ("human", "{input}")
]).partial(system_prompt=system_prompt)

llm = init_chat_model(model=MODEL_NAME, model_provider=MODEL_PROVIDER, base_url=BASE_URL, api_key=API_KEY, temperature=0.2)
chain = prompt | llm

msg = chain.invoke({"input": "What is 1 + 2?"})
raw = msg.content 
```

As we recall from the previous lecture, the model will return a string. We need to parse it to get the JSON object.

```python
# A registry is just a dictionary that lets us look up the right function by name
registry = {f.__name__: f for f in funcs}
call = json.loads(raw)
result = registry[call["tool"]](**call["args"]) if call["tool"] else None
print("Result:", result)
```

> [!NOTE]
> In practice, the model might return text around the JSON or malformed JSON. Robust implementations would add error handling here. We'll see in the exercises how to handle this.

If everything went well, we should have the result of the operation. In the next section we will see how to do this in a more general way with LangChain.

## Exercise

Time to put this into practice! See [exercise_manual_tools.md](exercises/0_exercise/exercise_manual_tools.md) — you will extend the manual calculator with a `divide` tool and add error handling.

