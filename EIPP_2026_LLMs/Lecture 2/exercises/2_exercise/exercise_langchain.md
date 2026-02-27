# Exercise 2 — LangChain Basics

## Goal

Practice building LangChain chains, working with different message types, and using prompt templates with `partial()` and `MessagesPlaceholder` — all inside your **`playground/`** environment.

## Prerequisites

- You have completed [Exercise 0](../0_exercise/exercise_uv.md) (playground setup)
- You have installed the LangChain packages (see [2_inference_with_langchain.md](../../2_inference_with_langchain.md))
- Your `.env` file contains API keys (configured in Exercise 1, though these exercises use the EBI cluster by default)

## Setup

### 1. Install the required packages

Make sure you are in the `playground/` directory and install the LangChain packages:

> [!TIP] Terminal
> ```bash
> uv add langchain langchain-google-genai langchain-openai python-dotenv
> ```

> [!NOTE]
> If you already ran this command while reading the guide, UV will skip packages that are already installed. This is normal.

### 2. Create the exercise folder and copy the starter scripts

> [!TIP] Terminal
> ```bash
> mkdir -p lecture_2/2_langchain
> cp "../Lecture 2/exercises/2_exercise/exercise_1_basic_chain.py" lecture_2/2_langchain/
> cp "../Lecture 2/exercises/2_exercise/exercise_2_conversation.py" lecture_2/2_langchain/
> cp "../Lecture 2/exercises/2_exercise/exercise_3_multi_specialist.py" lecture_2/2_langchain/
> ```

> [!NOTE]
> You can also copy the files using your operating system's file explorer instead of the terminal. Just drag or copy the three Python files from `Lecture 2/exercises/2_exercise/` into `playground/lecture_2/2_langchain/`.

All commands below assume you are in the `playground/` folder.

> [!NOTE]
> The exercise scripts are configured to use the EBI GPU cluster endpoint via ngrok. If the endpoint is unavailable, you can switch to Gemini by uncommenting the Gemini lines and commenting out the ngrok lines at the top of each script.

---

## Exercise 1: Basic Chain Construction (Easy)

Create a simple LangChain pipeline that takes a scientific concept and explains it in terms a high school student would understand.

**Your tasks** (edit `lecture_2/2_langchain/exercise_1_basic_chain.py`):

1. Create a `ChatPromptTemplate` using `from_messages` with:
   - A system message instructing the model to act as a friendly science teacher who explains complex topics simply
   - A human message with a `{concept}` placeholder

> [!TIP]
> Use `from_messages()` when you want separate system and human messages (like here). Use `from_template()` for simpler single-string prompts without a system message.
2. Build a chain using the pipe operator: `prompt | llm | StrOutputParser()`
3. Test it by invoking the chain with `{"concept": "photosynthesis"}`
4. Print the result

**Run your script:**

> [!TIP] Terminal
> ```bash
> uv run lecture_2/2_langchain/exercise_1_basic_chain.py
> ```

**Expected outcome**: A plain-text explanation of photosynthesis that a high school student could understand.

---

## Exercise 2: Message Types and Conversation (Easy)

Build a conversation that simulates a patient asking a doctor about medication side effects. This exercise practices using different message types and managing conversation history.

**Your tasks** (edit `lecture_2/2_langchain/exercise_2_conversation.py`):

1. Create a `SystemMessage` that makes the model act as a knowledgeable but empathetic doctor
2. Build a conversation with at least 2 exchanges:
   - Patient asks about side effects of a medication (e.g., ibuprofen)
   - Doctor responds (you invoke the model and capture the response)
   - Patient asks a follow-up question
   - Doctor responds again (using the full conversation history)
3. Print each response after receiving it
4. Print the type of each message in the conversation history at the end

**Run your script:**

> [!TIP] Terminal
> ```bash
> uv run lecture_2/2_langchain/exercise_2_conversation.py
> ```

**Expected outcome**: A realistic doctor-patient dialogue where the doctor remembers what was discussed earlier.

---

## Exercise 3: Multi-Specialist Consultation System (Hard)

Create a system where a patient can consult different medical specialists, each with their own expertise and conversation history.

**Your tasks** (edit `lecture_2/2_langchain/exercise_3_multi_specialist.py`):

1. Create a base `ChatPromptTemplate` with a `{specialty}` placeholder in the system message, a `MessagesPlaceholder` for history, and a human `{question}`
2. Use `partial()` on the base prompt to create prompts for 2 specialties (e.g., cardiology and neurology)
3. Create a `make_specialist_chain()` function that generates a chain for any given specialty

   **Hint:** Your function signature should look like:
   ```python
   def make_specialist_chain(specialty: str):
       # Use base_prompt.partial(specialty=specialty) | llm
       # Do NOT add StrOutputParser() — you need usage_metadata
       pass
   ```

4. Ask the same patient question to both specialists and print each response with token usage (`response.usage_metadata`)

> [!IMPORTANT]
> Do NOT use `StrOutputParser()` in your specialist chains if you need `usage_metadata` — the parser converts the response to a plain string, losing the metadata. Use `prompt | llm` without the parser, then access `response.usage_metadata` on the `AIMessage` object.

5. Pick one specialist, reuse its response from step 4 as conversation history, and ask a follow-up question to show the specialist remembers context

> [!NOTE]
> This exercise makes 3 API calls in total.

**Run your script:**

> [!TIP] Terminal
> ```bash
> uv run lecture_2/2_langchain/exercise_3_multi_specialist.py
> ```

**Expected outcome**: A system showing how different specialists respond to the same patient, with context maintained across exchanges.

---

## What your playground should look like now

After completing this exercise:

```
playground/
├── .env              ← your API key (never commit this!)
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
    └── 2_langchain/
        ├── exercise_1_basic_chain.py
        ├── exercise_2_conversation.py
        └── exercise_3_multi_specialist.py
```
