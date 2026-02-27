"""
Solution 0 — Manual Tool Calling
==================================
Extends the manual calculator with a divide tool and error handling.

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

# ── Tool definitions ─────────────────────────────────────

def add(a: float, b: float) -> float:
    """Add two numbers together"""
    return a + b

def multiply(a: float, b: float) -> float:
    """Multiply two numbers together"""
    return a * b

def divide(a: float, b: float):
    """Divide the first number by the second number"""
    if b == 0:
        return "Error: division by zero"
    return a / b

# ── Helper functions ─────────────────────────────────────

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

# --- Configuration ---
# MODEL_NAME = "gemini-3-flash-preview"
# MODEL_PROVIDER = "google_genai"
MODEL_NAME = "my-model"
MODEL_PROVIDER = "openai"
BASE_URL = "https://myung-courageous-gino.ngrok-free.dev/v1"
API_KEY = "not-needed"

funcs = [add, multiply, divide]

system_prompt = build_system_prompt(funcs)
registry = {f.__name__: f for f in funcs}

llm = init_chat_model(
    model=MODEL_NAME,
    model_provider=MODEL_PROVIDER,
    base_url=BASE_URL,
    api_key=API_KEY,
    temperature=0,
)

prompt = ChatPromptTemplate.from_messages([
    ("system", "{system_prompt}"),
    ("human", "{input}")
]).partial(system_prompt=system_prompt)

chain = prompt | llm

# ── Invoke and parse ─────────────────────────────────────

questions = [
    "What is 10 + 25?",
    "What is 7 * 8?",
    "What is 100 / 4?",
    "What is 5 / 0?",
]

for question in questions:
    print(f"\nQuestion: {question}")
    msg = chain.invoke({"input": question})
    raw = msg.content
    print(f"Raw model output: {raw}")

    try:
        call = json.loads(raw)
    except json.JSONDecodeError:
        print(f"Failed to parse JSON. Raw output was: {raw}")
        print("The model may have returned extra text around the JSON.")
        continue

    if call["tool"]:
        result = registry[call["tool"]](**call["args"])
        print(f"Tool: {call['tool']}")
        print(f"Args: {call['args']}")
        print(f"Result: {result}")
    else:
        print("No tool was called.")
