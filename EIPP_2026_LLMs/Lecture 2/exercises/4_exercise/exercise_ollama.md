# Exercise 4 — Ollama

## Goal

Practice using Ollama to run LLMs locally and integrate them with LangChain — all inside your **`playground/`** environment.

## Prerequisites

- You have completed the previous exercises
- You have installed Ollama and can run `ollama --version` in your terminal
- You have the `qwen3:0.6b` model available (either via `ollama pull qwen3:0.6b` or it was automatically downloaded when you ran `ollama run qwen3:0.6b` in the guide)
- You have read the guide [4_ollama.md](../../4_ollama.md)

## Setup

### 1. Install the required packages

Make sure you are in the `playground/` directory and install the Ollama integration for LangChain:

> [!TIP] Terminal
> ```bash
> uv add langchain-ollama
> ```

### 2. Create the exercise folder and copy the starter scripts

> [!TIP] Terminal
> ```bash
> mkdir -p lecture_2/4_ollama
> cp "../Lecture 2/exercises/4_exercise/check_ollama.py" lecture_2/4_ollama/
> cp "../Lecture 2/exercises/4_exercise/exercise_1_ollama_langchain.py" lecture_2/4_ollama/
> cp "../Lecture 2/exercises/4_exercise/exercise_2_local_vs_cloud.py" lecture_2/4_ollama/
> ```

> [!NOTE]
> You can also copy the files using your operating system's file explorer.

All commands below assume you are in the `playground/` folder.

### 3. Verify that Ollama is running

Before starting the exercises, run the verification script to make sure Ollama is reachable and the model is available:

> [!TIP] Terminal
> ```bash
> uv run lecture_2/4_ollama/check_ollama.py
> ```

If the check fails:
1. Make sure Ollama is installed: `ollama --version`
2. Make sure the service is running:
   - **macOS:** Open the Ollama app from Applications (or check if it's in the menu bar)
   - **Linux:** `sudo systemctl status ollama` — if not running, `sudo systemctl start ollama`
   - **Windows:** Check the system tray for the Ollama icon
3. Make sure the model is downloaded: `ollama list` should show `qwen3:0.6b`

---

## Exercise 1: Ollama with LangChain (Easy)

Use `init_chat_model` to connect LangChain to your local Ollama model and ask it a biology question.

**Your tasks** (edit `lecture_2/4_ollama/exercise_1_ollama_langchain.py`):

1. Create an LLM using `init_chat_model` with `model="qwen3:0.6b"` and `model_provider="ollama"`
2. Create a `ChatPromptTemplate` with a system message (act as a helpful biology tutor) and a human `{question}` placeholder
3. Build a chain: `prompt | llm | StrOutputParser()`
4. Invoke the chain with a biology question of your choice
5. Print the result

**Run your script:**

> [!TIP] Terminal
> ```bash
> uv run lecture_2/4_ollama/exercise_1_ollama_langchain.py
> ```

**Expected outcome**: A text explanation of your biology topic, generated entirely on your local machine.

> [!NOTE]
> The `qwen3:0.6b` model is very small (0.6 billion parameters), so answers will be **shorter and less detailed** than what you'd get from Gemini or larger cloud models. This is expected — the trade-off for running locally in under 600 MB of memory.

> [!NOTE]
> The first time you invoke a model through LangChain, it may take **10–30 seconds** to load into memory. You'll see no output during this time — **don't cancel the script!** Subsequent calls will be much faster (1–2 seconds).

---

## Exercise 2: Local vs Cloud Comparison (Medium)

Compare a locally-deployed Ollama model with a cloud-based model on the same question. Measure response time and compare output quality.

**Your tasks** (edit `lecture_2/4_ollama/exercise_2_local_vs_cloud.py`):

1. Create a local LLM using `init_chat_model` with `model_provider="ollama"` and `model="qwen3:0.6b"`
2. Create a cloud LLM using `init_chat_model` with `model_provider=MODEL_PROVIDER` and `model=MODEL_NAME`

> [!NOTE]
> The starter script defines `MODEL_PROVIDER` and `MODEL_NAME` at the top — this is the same configuration block you've seen in previous exercises. Use these variables when creating the cloud LLM.
3. Define a shared question (e.g., "Explain how vaccines work in simple terms.")
4. For each model:
   - Record the start time with `time.time()`
   - Invoke the model with the question
   - Record the end time
   - Print the model name, response text, and elapsed time
5. Print a separator between models for easy comparison

**Run your script:**

> [!TIP] Terminal
> ```bash
> uv run lecture_2/4_ollama/exercise_2_local_vs_cloud.py
> ```

**Expected outcome**: Side-by-side comparison of a local model and a cloud model, with timing data showing the trade-offs.

> [!NOTE]
> The local model (qwen3:0.6b) is very small and will produce shorter, less detailed answers than the cloud model. This is expected — the point of this exercise is to see the practical differences between local and cloud inference.

> [!NOTE]
> The cloud part of this exercise uses the EBI GPU cluster endpoint via ngrok. If the endpoint is unavailable, you can switch to Gemini by uncommenting the Gemini lines and commenting out the ngrok lines at the top of the script.

---

## What your playground should look like now

After completing this exercise:

```
playground/
├── .env              ← your API keys (never commit this!)
├── .gitignore
├── .venv/
├── pyproject.toml
├── uv.lock
└── lecture_2/
    ├── 0_uv_setup/
    │   └── check_uv_setup.py
    ├── 1_gemini_api/
    │   ├── check_api_key.py
    │   └── gemini_biology.py
    ├── 2_langchain/
    │   ├── exercise_1_basic_chain.py
    │   ├── exercise_2_conversation.py
    │   └── exercise_3_multi_specialist.py
    ├── 3_openrouter/
    │   ├── exercise_1_openrouter_call.py
    │   └── exercise_2_model_comparison.py
    └── 4_ollama/
        ├── check_ollama.py
        ├── exercise_1_ollama_langchain.py
        └── exercise_2_local_vs_cloud.py
```
