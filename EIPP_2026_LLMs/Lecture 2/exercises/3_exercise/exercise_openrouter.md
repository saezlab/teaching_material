# Exercise 3 — OpenRouter

## Goal

Practice using OpenRouter to access different LLM providers through a single API, compare models, and integrate OpenRouter models with LangChain — all inside your **`playground/`** environment.

## Prerequisites

- You have completed [Exercise 0](../0_exercise/exercise_uv.md), [Exercise 1](../1_exercise/exercise_gemini_api.md), and [Exercise 2](../2_exercise/exercise_langchain.md)
- Your `.env` file in `playground/` contains your `OPENROUTER_API_KEY`
- You have read the guide [3_extension_to_openrouter.md](../../3_extension_to_openrouter.md)

## Step 0: Verify Your Setup

Before starting, confirm your OpenRouter API key is loaded correctly:

> [!TIP] Terminal
> ```bash
> mkdir -p lecture_2/3_openrouter
> cp "../Lecture 2/exercises/3_exercise/check_openrouter_api_key.py" lecture_2/3_openrouter/
> uv run lecture_2/3_openrouter/check_openrouter_api_key.py
> ```

If this fails, check your `.env` file in `playground/`.

## Setup

### 1. Install the required packages

Make sure you are in the `playground/` directory and install the OpenAI-compatible LangChain package:

> [!TIP] Terminal
> ```bash
> uv add langchain-openai
> ```

### 2. Create the exercise folder and copy the starter scripts

> [!TIP] Terminal
> ```bash
> mkdir -p lecture_2/3_openrouter
> cp "../Lecture 2/exercises/3_exercise/exercise_1_openrouter_call.py" lecture_2/3_openrouter/
> cp "../Lecture 2/exercises/3_exercise/exercise_2_model_comparison.py" lecture_2/3_openrouter/
> ```

> [!NOTE]
> You can also copy the files using your operating system's file explorer instead of the terminal. Just drag or copy the two Python files from `Lecture 2/exercises/3_exercise/` into `playground/lecture_2/3_openrouter/`.

All commands below assume you are in the `playground/` folder.

---

## Exercise 1: First OpenRouter Call (Easy)

Create a model that connects to OpenRouter using `init_chat_model`, and ask it a biology question.

**Your tasks** (edit `lecture_2/3_openrouter/exercise_1_openrouter_call.py`):

1. Create an LLM using `init_chat_model` with:
   - `model="google/gemma-3-27b-it:free"`
   - `model_provider="openai"`
   - `base_url="https://openrouter.ai/api/v1"`
   - `api_key=os.environ["OPENROUTER_API_KEY"]`
2. Invoke the model with a biology question (e.g., "What are the main differences between DNA and RNA?")
3. Print the response text with `response.text`
4. Print the token usage with `response.usage_metadata`

**Run your script:**

> [!TIP] Terminal
> ```bash
> uv run lecture_2/3_openrouter/exercise_1_openrouter_call.py
> ```

**Expected outcome**: A text response from the model + token usage stats.

> [!NOTE]
> Some free models on OpenRouter report `output_tokens: 0` in the usage metadata. This is a known quirk of certain free-tier providers — your response still arrived, the token count just isn't tracked.

---

## Exercise 2: Model Comparison (Medium)

Compare two free OpenRouter models on the same biology question. Measure quality and speed.

**Your tasks** (edit `lecture_2/3_openrouter/exercise_2_model_comparison.py`):

1. Create a helper function `make_openrouter_llm(model_name, temperature=0.2)` that returns an `init_chat_model` instance configured for OpenRouter
2. Define a list of 2 free models (e.g., `["google/gemma-3-27b-it:free", "deepseek/deepseek-r1-0528:free"]`)
3. For each model:
   - Measure response time with `time.time()`
   - Invoke the model with the same question
   - Print the model name, response text, token usage, and elapsed time
4. Print a separator line between models so you can easily compare

> [!TIP]
> **Suggested model pairs for comparison:**
> - Similar models: `["google/gemma-3-27b-it:free", "google/gemma-3-9b-it:free"]`
> - Different families: `["google/gemma-3-27b-it:free", "meta-llama/llama-3.2-3b-instruct:free"]`

> [!NOTE]
> Free models on OpenRouter change over time. If a model in your list is unavailable, browse https://openrouter.ai/models?pricing=free to find alternatives.

> [!NOTE]
> Reasoning models (like `deepseek/deepseek-r1-0528:free`) may take significantly longer to respond because they "think" before answering. You may also see a `reasoning` field in the token usage — those are internal thinking tokens. This is normal and makes for an interesting comparison!

> [!TIP]
> If you see a `404` error or "model not found," the model may have been removed from OpenRouter's free tier. Browse https://openrouter.ai/models?pricing=free to find current alternatives.

**Run your script:**

> [!TIP] Terminal
> ```bash
> uv run lecture_2/3_openrouter/exercise_2_model_comparison.py
> ```

**Expected outcome**: Side-by-side comparison of model responses with timing data.

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
    └── 3_openrouter/
        ├── check_openrouter_api_key.py
        ├── exercise_1_openrouter_call.py
        └── exercise_2_model_comparison.py
```
