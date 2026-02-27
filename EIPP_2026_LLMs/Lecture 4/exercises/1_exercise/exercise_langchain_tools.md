# Exercise 1 — LangChain Tools

## Goal

Practice defining tools with LangChain's `@tool` decorator, binding them to an LLM, executing tool calls, and returning results for a final answer.

## Prerequisites

- You have completed Lecture 2 exercises and Lecture 3 exercises

## Setup

### 1. Create the exercise folder and copy the starter script

> [!TIP] Terminal
> ```bash
> mkdir -p lecture_4/1_langchain_tools
> cp "../Lecture 4/exercises/1_exercise/exercise_langchain_tools.py" lecture_4/1_langchain_tools/
> ```

> [!NOTE]
> You can also copy the file using your file explorer. Just drag the Python file from `Lecture 4/exercises/1_exercise/` into `playground/lecture_4/1_langchain_tools/`.

All commands below assume you are in the `playground/` folder.

---

## The Exercise

In this exercise you will define two biology-themed tools and have an LLM use them to answer questions about DNA sequences.

**Your tasks** (edit `lecture_4/1_langchain_tools/exercise_langchain_tools.py`):

1. **Define a `melting_temperature` tool** that estimates the melting temperature of a short DNA sequence using the Wallace rule. The Wallace rule is a simple approximation for estimating the melting temperature of short DNA oligonucleotides (<=14 bases):
   - Tm = 2 * (A+T count) + 4 * (G+C count)
   - For this exercise, just apply the formula regardless of length. In practice, more accurate methods (like nearest-neighbor) are used for longer sequences.
   - The function should take a `sequence: str` and return a `float`
   - Add the `@tool` decorator and a clear docstring

2. **Define a `reverse_complement` tool** that returns the reverse complement of a DNA sequence:
   - A↔T, C↔G
   - The function should take a `sequence: str` and return a `str`
   - Add the `@tool` decorator and a clear docstring

3. **Bind both tools** to the LLM using `llm.bind_tools([...])`

4. **Execute the tool calls**: When the model returns tool calls, look up each tool by name and invoke it with the provided arguments

5. **Return results to the model**: Send the tool results back as `ToolMessage` objects so the model can write a final answer. Use the base `llm` (not `llm_with_tools`) for the final invocation so the model writes a text answer instead of potentially calling more tools.

> [!TIP] Terminal
> ```bash
> uv run lecture_4/1_langchain_tools/exercise_langchain_tools.py
> ```

**Expected outcome:** The model calls the appropriate tool(s) and then produces a final answer that includes the computed values (e.g., "The melting temperature of ATGCGATCG is approximately 28 degrees Celsius").

---

## Bonus Challenge

Try asking a question that requires **both** tools, for example:
> "What is the reverse complement of ATGCGATCG and what is the melting temperature of the original sequence?"

Does the model call both tools? How does it combine the results?

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
│   └── 1_langchain_tools/
│       └── exercise_langchain_tools.py
```
