"""
Exercise 0 — Manual Tool Calling
=================================
Extend the manual calculator with a divide tool and add error handling
for JSON parsing.

Your tasks:
  1. Add a divide(a, b) function that handles division by zero
  2. Register it in the funcs list
  3. Add try/except around JSON parsing
  4. Test with different questions

Usage (from the playground/ root):
    uv run lecture_4/0_manual_tools/exercise_manual_tools.py
"""

from dotenv import load_dotenv
load_dotenv()

import json
import inspect
from typing import Callable, List, get_type_hints

from langchain.chat_models import init_chat_model
from langchain_core.prompts import ChatPromptTemplate

# --- Configuration ---
# MODEL_NAME = "gemini-3-flash-preview"
# MODEL_PROVIDER = "google_genai"
MODEL_NAME = "my-model"
MODEL_PROVIDER = "openai"
BASE_URL = "https://myung-courageous-gino.ngrok-free.dev/v1"
API_KEY = "not-needed"

# ── Tool definitions ─────────────────────────────────────

def add(a: float, b: float) -> float:
    """Add two numbers together"""
    return a + b

def multiply(a: float, b: float) -> float:
    """Multiply two numbers together"""
    return a * b

# TODO: Define a divide(a: float, b: float) function here.
# It should return a / b, but handle division by zero gracefully
# (return an error string like "Error: division by zero").
# Don't forget the docstring!


# ── Helper functions (already provided — no changes needed) ───

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

# Simplified version of the guide's build_system_prompt
def build_system_prompt(funcs: List[Callable]) -> str:
    tools_block = build_tools_block(funcs)
    schema = build_json_schema(funcs)
    return f"""You can call Python tools by emitting exactly ONE JSON object.

- tool: one of the tool names below (or null if no tool is needed)
- args: an object whose keys match the tool's parameters exactly

TOOLS AVAILABLE:
{tools_block}

OUTPUT FORMAT (STRICT):
Return ONLY a JSON object, no prose, no code fences.

Call object schema:
{schema}

If you do NOT need any tool, return:
{{"tool": null, "args": {{}}}}
"""

# ── Build the chain ──────────────────────────────────────

# TODO: Add your divide function to this list
funcs = [add, multiply]

system_prompt = build_system_prompt(funcs)
registry = {f.__name__: f for f in funcs}

llm = init_chat_model(
    model=MODEL_NAME,
    model_provider=MODEL_PROVIDER,
    base_url=BASE_URL,
    api_key=API_KEY,
    temperature=0.2,
)

prompt = ChatPromptTemplate.from_messages([
    ("system", "{system_prompt}"),
    ("human", "{input}")
]).partial(system_prompt=system_prompt)

chain = prompt | llm

# ── Invoke and parse ─────────────────────────────────────

question = "What is 10 + 25?"
print(f"Question: {question}")

msg = chain.invoke({"input": question})
raw = msg.content
print(f"Raw model output: {raw}")

# TODO: Wrap this section in a try/except to catch json.JSONDecodeError.
# If parsing fails, print a helpful message instead of crashing.
call = json.loads(raw)

if call["tool"]:
    result = registry[call["tool"]](**call["args"])
    print(f"Tool: {call['tool']}")
    print(f"Args: {call['args']}")
    print(f"Result: {result}")
else:
    print("No tool was called.")
