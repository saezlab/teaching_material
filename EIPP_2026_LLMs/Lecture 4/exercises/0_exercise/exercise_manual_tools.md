# Exercise 0 — Manual Tool Calling

## Goal

Practice the from-scratch approach to tool calling: define Python functions, describe them to the model via a system prompt, parse the model's JSON output, and execute the function.

## Prerequisites

- You have read [`0_introduction_to_tools.md`](../../0_introduction_to_tools.md)
- You have completed Lecture 2 exercises (playground setup, Gemini API, LangChain basics) and Lecture 3 exercises

## Setup

### 1. Create the exercise folder and copy the starter script

> [!TIP] Terminal
> ```bash
> mkdir -p lecture_4/0_manual_tools
> cp "../Lecture 4/exercises/0_exercise/exercise_manual_tools.py" lecture_4/0_manual_tools/
> ```

> [!NOTE]
> You can also copy the file using your file explorer. Just drag the Python file from `Lecture 4/exercises/0_exercise/` into `playground/lecture_4/0_manual_tools/`.

All commands below assume you are in the `playground/` folder.

---

## The Exercise

The starter script already provides:
- `add(a, b)` and `multiply(a, b)` tool functions
- Helper functions that build a system prompt with tool descriptions and JSON schema
- A working example that invokes the model and parses the JSON response

**Your tasks** (edit `lecture_4/0_manual_tools/exercise_manual_tools.py`):

1. **Add a `divide` tool** that divides two numbers. It should:
   - Take two `float` arguments (`a` and `b`)
   - Return `a / b`
   - Handle division by zero by returning an error string like `"Error: division by zero"`
   - Have a clear docstring

2. **Register the new tool** in the `funcs` list so the system prompt includes it

3. **Add error handling** to the JSON parsing step:
   - Wrap `json.loads(raw)` in a try/except that catches `json.JSONDecodeError`
   - If parsing fails, print the raw output and a helpful message

   Hint: Your except block might look like:
   ```python
   except json.JSONDecodeError:
       print(f"Failed to parse JSON. Raw output was: {raw}")
   ```

4. **Test your tools** by changing the input question to test:
   - An addition (e.g., "What is 10 + 25?")
   - A multiplication (e.g., "What is 7 * 8?")
   - A division (e.g., "What is 100 / 4?")
   - A division by zero (e.g., "What is 5 / 0?")

> [!TIP] Terminal
> ```bash
> uv run lecture_4/0_manual_tools/exercise_manual_tools.py
> ```

> [!TIP] LLM Smart Behavior
> When you test "5 / 0", the LLM might be smart enough to recognize division by zero and **choose not to call the divide tool at all**, returning `{"tool": null, "args": {}}` instead. This is actually interesting behavior — the model is avoiding an error proactively. If this happens, you can still see it's working correctly.

**Expected outcome:** The script correctly calls the right tool for each question. Division by zero should return the error string instead of crashing.

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
└── lecture_4/
    └── 0_manual_tools/
        └── exercise_manual_tools.py
```
