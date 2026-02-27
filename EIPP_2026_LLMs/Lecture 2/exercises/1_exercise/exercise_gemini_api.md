# Exercise: Using the Gemini API

## Goal

Set your Gemini API key as an environment variable, verify it works, and run a Python script that sends a biological text to Gemini for analysis — all inside your **`playground/`** environment.

## Prerequisites

- You have completed [Exercise 0](../0_exercise/exercise_uv.md) and your `playground/` is initialized with UV
- You have obtained a Gemini API key (see [1_getting_gemini_API_key.md](../../1_getting_gemini_API_key.md))

## Steps

### 1. Navigate to the playground and create the exercise folder

Make sure you are in the `playground/` directory and create a subfolder for this exercise:

> [!TIP] Terminal
> ```bash
> $ cd playground
> $ mkdir -p lecture_2/1_gemini_api
> ```

> [!NOTE]
> You can also create the folders using your operating system's file explorer (Finder on macOS, File Explorer on Windows) instead of the terminal.

### 2. Copy the exercise scripts

> [!TIP] Terminal
> ```bash
> $ cp ../Lecture\ 2/exercises/1_exercise/check_api_key.py lecture_2/1_gemini_api/
> $ cp ../Lecture\ 2/exercises/1_exercise/gemini_biology.py lecture_2/1_gemini_api/
> ```

> [!NOTE]
> You can also copy the files using your operating system's file explorer instead of the terminal. Just drag or copy `check_api_key.py` and `gemini_biology.py` from `Lecture 2/exercises/1_exercise/` into `playground/lecture_2/1_gemini_api/`.

### 3. Set your API key as an environment variable

The standard way to pass secrets (like API keys) to a program is through **environment variables**. This avoids hardcoding keys in your source code, which is a security risk — especially if you push your code to a public repository.

Set the variable in your current terminal session:

**macOS/Linux**:

> [!TIP] Terminal
> ```bash
> export GOOGLE_API_KEY='your-api-key-here'
> ```

**Windows (PowerShell)**:

> [!TIP] Terminal
> ```powershell
> $env:GOOGLE_API_KEY = "your-api-key-here"
> ```

> **Note**: This only lasts for the current terminal session. If you close the terminal, you will need to set it again.

### 4. Create a `.env` file to persist your API key

Setting the environment variable manually every time you open a new terminal is tedious. A common solution is to store the key in a **`.env` file** at the root of your project. Future exercises will load it automatically using the `python-dotenv` package.

Create the file in your `playground/` directory:

> [!TIP] Terminal
> ```bash
> $ echo 'GOOGLE_API_KEY=your-api-key-here' > .env
> ```

Replace `your-api-key-here` with your actual API key (no quotes around the value needed in the file).

> [!NOTE]
> You can also create the `.env` file manually: open a text editor, type `GOOGLE_API_KEY=your-api-key-here`, and save it as `.env` in the `playground/` folder.

Your `.env` file should look like this:

```
GOOGLE_API_KEY=AIzaSy...your-key-here
```

#### Keep your `.env` file out of version control

The `.env` file contains a secret. If you ever push your `playground/` to GitHub, you must make sure it is **not** included. Add it to `.gitignore`:

> [!TIP] Terminal
> ```bash
> $ echo '.env' >> .gitignore
> ```

> [!IMPORTANT]
> Never commit your `.env` file to a git repository. If you accidentally push an API key to a public repository, anyone can use (and abuse) it.

#### How future exercises will use this file

In subsequent exercises, scripts will load the `.env` file automatically using `python-dotenv`. The typical pattern looks like this:

```python
from dotenv import load_dotenv
load_dotenv()  # reads .env and sets the variables in os.environ

import os
api_key = os.environ.get("GOOGLE_API_KEY")
```

You don't need to worry about this yet — the current exercise uses the environment variable you exported in Step 3. But from the next exercise onward, scripts will pick up the key from `.env` so you won't need to `export` it manually anymore.

> [!IMPORTANT]
> The verification script reads the API key from your terminal session (via `export`), **not** from the `.env` file. Make sure you ran the `export` command from Step 3 in the same terminal where you run the script below. The `.env` file will be used automatically from Exercise 2 onward.

### 5. Verify your API key

Run the verification script from the `playground/` root:

> [!TIP] Terminal
> ```bash
> $ uv run lecture_2/1_gemini_api/check_api_key.py
> ```

If everything is set up correctly, you should see Gemini respond with a confirmation message. If it fails, double-check that:

- You copied the full API key without extra spaces
- You exported it in the **same** terminal where you run the script

### 6. Install the `google-genai` package

The biology script needs the `google-genai` package to talk to the Gemini API. Add it to your project just like you did with `cowsay` in Exercise 0:

> [!TIP] Terminal
> ```bash
> $ uv add google-genai
> ```

This permanently adds the package to your `pyproject.toml` so it is available for all future exercises that use the Gemini API.

### 7. Run the biology exercise

Now you can run the script with a simple `uv run`:

> [!TIP] Terminal
> ```bash
> $ uv run lecture_2/1_gemini_api/gemini_biology.py
> ```

The script sends a paragraph about the CRISPR-Cas9 genome editing system to Gemini and asks it to:

1. Extract the key biological concepts
2. Explain each one for an undergraduate audience
3. Suggest a follow-up question

> [!NOTE]
> **Rate limiting:** Since we are using free models, you may occasionally get rate-limited (you will see an error containing `429 RESOURCE_EXHAUSTED`). The traceback can look long and intimidating, but the important part is at the bottom. If this happens, wait about 60 seconds and try again. If the problem persists, you can change the model from `gemini-3-flash-preview` to `gemini-2.5-flash` in the script — this uses a different quota and should work immediately.

### 8. Experiment

Try modifying `lecture_2/1_gemini_api/gemini_biology.py`:

- Change the `BIOLOGICAL_TEXT` to a different topic (e.g. protein folding, PCR, gene expression)
- Change the `PROMPT` to ask Gemini something different (e.g. "translate this to plain English", "generate quiz questions from this text")
- Try a different model by changing the value of `MODEL_NAME` to another available model

## How environment variables work

When you run `export GOOGLE_API_KEY='...'`, the key is stored in your shell's environment. Any program you launch from that shell can read it — in Python via `os.environ.get("GOOGLE_API_KEY")`.

This is the standard approach for managing secrets because:

- The key never appears in your source code
- It is not committed to version control (unlike a `.env` file, which can be accidentally pushed)
- It works the same way in local development, CI/CD, and production servers

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
    └── 1_gemini_api/
        ├── check_api_key.py
        └── gemini_biology.py
```
