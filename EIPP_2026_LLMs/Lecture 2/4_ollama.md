# Deploying an LLM locally with Ollama

## Table of Contents

1. [Deploying an LLM locally with Ollama](#deploying-an-llm-locally-with-ollama)
2. [Installing Ollama](#installing-ollama)
3. [Running models from the CLI](#running-models-from-the-cli)
4. [Using Ollama with LangChain](#using-ollama-with-langchain)
5. [Exercises](#exercises)

## Deploying an LLM locally with Ollama

**Ollama** is an **open-source tool** that lets you **run large language models (LLMs) locally** on your own computer—whether you're using Windows, macOS, or Linux—without needing to rely on cloud-based services.

Here are some of the main features of Ollama:

* **Privacy & Control:** Since everything runs on your device, your data stays local—offering stronger privacy compared to cloud-hosted models.
* **Efficiency:** Ollama uses techniques like **quantization** to minimize resource consumption, making LLMs viable on consumer-grade hardware.
* **Simplicity:** It provides a **command-line interface (CLI)** (and recently, a **Windows GUI**) that simplifies tasks like pulling, running, and managing models.
* **Cross-Platform & Extensible:** Works across major operating systems and supports integration via APIs and model customizations (like defining prompts or using LoRA fine-tuning).

> [!WARNING]
> **Hardware Requirements**
>
> Local LLM deployment typically requires an NVIDIA GPU with 4–80GB+ VRAM depending on model size. However, **for this tutorial, we use the very small `qwen3:0.6b` model (523 MB), which runs fine on CPU** for basic experimentation.
>
> **Expected performance on CPU:** 10–30 seconds for the first response (model loading), 2–5 seconds for subsequent calls. This is slower than GPU (<1 second) but perfectly usable for learning.
>
> Later in the course, we'll explore running larger models on dedicated HPC infrastructure with appropriate GPU resources.

## Installing Ollama

**On Linux**, run the following command in your terminal:

```bash
curl -fsSL https://ollama.com/install.sh | sh
```

**On macOS**, you can install Ollama using Homebrew:

```bash
brew install ollama
```

Alternatively, you can download the official installer from [https://ollama.com/download](https://ollama.com/download).

**On Windows**, download the official installer from [https://ollama.com/download](https://ollama.com/download) and run it. After installation, the Ollama service starts automatically.

**Verify installation:**

> [!TIP] Terminal
> ```bash
> ollama --version
> # Should print: ollama version 0.x.x
> ```

**Check if the service is running:**

> [!TIP] Terminal
> ```bash
> curl http://localhost:11434
> # Should print: "Ollama is running"
> ```

If the curl command fails, start the Ollama service:
- **macOS:** Open the Ollama app from Applications (or check if it's in the menu bar)
- **Linux:** The service usually auto-starts; if not, run `sudo systemctl start ollama`
- **Windows:** The Ollama app should be running in the system tray

## Running models from the CLI

Once Ollama is installed and running, we are ready to deploy an LLM. A list of available models can be found [here](https://ollama.com/library).

Model names on Ollama follow two patterns:
- **`namespace/model`** for community-contributed models (e.g., `mapler/gpt2` — `mapler` is the contributor)
- **`model:tag`** for official models (e.g., `qwen3:0.6b` where `0.6b` = 0.6 billion parameters)

Let's start with GPT-2 (2019), one of the first widely-available generative language models. At 328 MB it's tiny by modern standards, and since it predates instruction fine-tuning, its responses will feel very different from today's conversational models. This makes it a great reference point for seeing how far the field has come:

```bash
ollama run mapler/gpt2
```

See if you can get the model to explain to you what the key messages of the 2017 paper "Attention is all you need" are. Or how it responds to you saying "Hi!".

> [!TIP]
> For understanding the difference between legacy models like GPT-2 and the
> conversational style of current models, we can look into advances in
> *instruction fine-tuning* and *reinforcement learning from human feedback*;
> see for example this [IBM blog post](https://www.ibm.com/topics/instruction-tuning).

We can now compare it with the latest version of the small qwen3 model. This 0.6 billion parameter model only has a size of 523 MB, and interacting with it should be fast. Also this model is quite different from SOTA commercial models available today but should deliver different results compared to GPT-2. Give it a try:

```bash
ollama run qwen3:0.6b
```

> [!NOTE]
> The first time you run a model, Ollama downloads it automatically (523 MB for `qwen3:0.6b`). This may take a few minutes depending on your internet connection. Subsequent runs will start instantly.

I encourage you to try the same queries as with the GPT-2 model and see how the responses differ.

When you're done chatting, type `/bye` to exit the chat session. Commands starting with `/` are Ollama-specific shortcuts — type `/help` to see all available commands.

## Using Ollama with LangChain

The LangChain library provides a wrapper for Ollama that allows you to use Ollama as the backend for LLM calls. First, make sure you are in the `playground/` folder and install the Ollama integration for LangChain:

> [!TIP] Terminal
> ```bash
> uv add langchain-ollama
> ```

Once you have the package installed and an Ollama model running, you can use `init_chat_model` with the proper arguments to tell LangChain how to connect to the model:

> [!NOTE]
> Unlike Gemini and OpenRouter, Ollama runs entirely on your local machine — **no API key or `load_dotenv()` needed**. The model communicates via `localhost`.

```python
from langchain.chat_models import init_chat_model

llm = init_chat_model(
    model="qwen3:0.6b",
    model_provider="ollama",
    temperature=0.2,
)

response = llm.invoke("What is the capital of France?")
print(response.text)
```

## Exercises

The exercises for this section are in the `exercises/4_exercise/` folder. See [exercise_ollama.md](exercises/4_exercise/exercise_ollama.md) for instructions.

The exercises will guide you through:
1. **Exercise 1** (Easy): Use Ollama with LangChain to query a local model
2. **Exercise 2** (Medium): Compare a local model with a cloud model
