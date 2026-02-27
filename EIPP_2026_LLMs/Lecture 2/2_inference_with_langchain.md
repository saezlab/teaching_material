# Index

* [Installing the required packages](#installing-the-required-packages)
* [What is LangChain?](#what-is-langchain)
* [One-liner init + first invoke](#one-liner-init--first-invoke)
* [Different types of messages](#different-types-of-messages)
* [Statelessness](#statelessness)
* [Dissecting the structure of the response](#dissecting-the-structure-of-the-response)
* [The pipe operator in LangChain](#the-pipe-operator-in-langchain)
* [Reusing prompts with partial() syntax](#reusing-prompts-with-partial-syntax)
* [Conversation memory via `MessagesPlaceholder`](#conversation-memory-via-messagesplaceholder)
* [Swapping providers without refactoring the chains](#swapping-providers-without-refactoring-the-chains)
* [Exercises](#exercises)

## Installing the required packages

Before we start, we need to install the LangChain packages into our playground environment. Make sure you are in the `playground/` folder, then run:

> [!TIP] Terminal
> ```bash
> uv add langchain langchain-google-genai langchain-openai python-dotenv
> ```

This installs four packages:
- **langchain** — the core framework (also pulls in `langchain-core` automatically)
- **langchain-google-genai** — connects LangChain to the Gemini API
- **langchain-openai** — connects LangChain to OpenAI-compatible APIs (used for the EBI GPU cluster endpoint)
- **python-dotenv** — lets you load your API key from the `.env` file
After running the command, check your `pyproject.toml` to confirm the new dependencies appear.

## What is LangChain?

LangChain is a **framework for building applications with large language models (LLMs)**. Instead of treating an LLM as a black box that takes a prompt and returns text, LangChain provides a structured way to **connect models with other components** like prompts, memory, external APIs, databases, and tools.

The key ideas are:

* **Unified interface across LLM providers** (companies like Google, OpenAI, Meta that offer model APIs): You can call models from Google (Gemini), OpenAI, Anthropic, etc., with the same syntax, making it easier to switch providers.
* **Messages and prompts as building blocks**: Instead of writing plain strings, LangChain uses structured message objects (`SystemMessage`, `HumanMessage`, `AIMessage`) and templates (`ChatPromptTemplate`) to standardize conversations.
* **Chains and runnables**: You can link multiple steps—like prompt construction → model call → output parsing—into a **chain**, making LLM workflows more reusable and testable.
* **Memory and placeholders**: LangChain supports injecting **chat history** and **variables** dynamically into prompts, enabling contextual conversations.
* **Extensibility**: Beyond calling models, LangChain is designed to integrate with search engines, knowledge bases, and other tools, letting you build **reasoning-and-action pipelines** like RAG (**R**etrieval-**A**ugmented **G**eneration — connecting models to external knowledge bases) or agent systems.

In short:
LangChain is not just about *talking* to a model—it's about **orchestrating LLMs together with prompts, context, and external resources** to build reliable applications.

For more information, see the [official documentation](https://python.langchain.com/docs/introduction/).


## One-liner init + first invoke

Let's start with a simple example of how to use LangChain to call a model. The `init_chat_model` function is used to initialize a chat model. The model is picked via the `model` and `model_provider` parameters. The `temperature` parameter controls the randomness of the model's output (`0.0` = deterministic/repeatable, `1.0` = creative/varied). For most technical tasks like biology Q&A, `0.2` is a good default.

> [!WARNING]
> **Loading Environment Variables**
>
> Before running any LangChain code, make sure your API keys are loaded. Add these lines at the top of every script:
> ```python
> from dotenv import load_dotenv
> load_dotenv()
> ```
> This loads your `.env` file from the current directory. You can verify the location with `import os; print(os.getcwd())`.

> [!NOTE]
> **About the Configuration Block**
>
> You'll notice two configurations in the examples below:
> - **EBI GPU cluster** (active by default): Uses a dedicated GPU server hosted at EBI via an ngrok tunnel. This provides faster, more reliable access without rate limits.
> - **Gemini API** (commented out): The fallback option using your Google API key.
>
> **Why use the EBI cluster?** Gemini's free tier has strict rate limits (15 requests/minute). For this course, we use the EBI cluster to avoid hitting these limits during exercises.
>
> **Fallback:** If the EBI endpoint is unavailable (e.g., maintenance), uncomment the Gemini lines and comment out the EBI lines.

```python
from langchain.chat_models import init_chat_model

# --- Configuration ---
# MODEL_NAME = "gemini-3-flash-preview"
# MODEL_PROVIDER = "google_genai"
MODEL_NAME = "my-model"
MODEL_PROVIDER = "openai"
BASE_URL = "https://myung-courageous-gino.ngrok-free.dev/v1"
API_KEY = "not-needed"

# Pick via provider:model string → works across providers
llm = init_chat_model(model=MODEL_NAME,
                      model_provider=MODEL_PROVIDER,
                      base_url=BASE_URL,
                      api_key=API_KEY,
                      temperature=0.2)

response = llm.invoke("What is the capital of France?")
print(response.text)
```

By using the `invoke` method, we can call the model with a prompt. The model will return a response. `init_chat_model` is a convenience function that allows to initialize different models from different providers with a unified interface. To see the list of supported models, see the [official documentation](https://python.langchain.com/api_reference/langchain/chat_models/langchain.chat_models.base.init_chat_model.html).


## Different types of messages

### Why messages instead of plain strings?

When you use ChatGPT or Gemini in the browser, the conversation looks like a simple back-and-forth. But behind the scenes, every API call sends a **structured list of messages**, each tagged with a **role**. This is how the model knows who said what:

| Role | LangChain class | Purpose |
|------|----------------|---------|
| **system** | `SystemMessage` | Sets the model's persona and rules. The model treats this as background instructions, not as part of the conversation. |
| **human** (or "user") | `HumanMessage` | What the user says — your questions, prompts, or input data. |
| **ai** (or "assistant") | `AIMessage` | The model's replies. When building conversation history, you include past AI responses so the model has context. |

A typical API call looks like this conceptually:

```
[
  SystemMessage: "You are a helpful biology tutor."     ← sets behavior
  HumanMessage:  "What is apoptosis?"                   ← user question
  AIMessage:     "Apoptosis is programmed cell death..." ← previous model reply
  HumanMessage:  "How does it differ from necrosis?"    ← new user question
]
```

The model reads this list top to bottom and generates the next `AIMessage`. **The roles matter** — during training, the model learned to behave differently based on whether text comes from the system, the user, or a previous assistant turn. If you mislabel messages, the model can give unexpected results.

> [!NOTE]
> LangChain handles the provider-specific formatting for you. Different models (Gemini, LLaMA, GPT, etc.) use different internal templates to separate roles, but you always use the same `SystemMessage` / `HumanMessage` / `AIMessage` classes.

**When to use each type:**
- Use **SystemMessage** to set persistent behavior across all interactions (e.g., "always cite sources," "explain in simple terms")
- Use **HumanMessage** for specific questions or inputs
- Include **AIMessage** when building conversation history to give the model context of previous exchanges

### Using message types in practice

```python
from langchain_core.messages import SystemMessage, HumanMessage, AIMessage

# helpful assistant
system_message = SystemMessage(content="You are a helpful assistant.")
human_message = HumanMessage(content="What is the capital of France?")
resp = llm.invoke([system_message, human_message])
print(resp.text)

# pirate assistant — same question, different system prompt
system_message = SystemMessage(content="You always answer like a pirate")
human_message = HumanMessage(content="What is the capital of France?")
resp = llm.invoke([system_message, human_message])
print(resp.text)

# check the type of all the messages
print(type(system_message))
print(type(human_message))
print(type(resp))
```

Notice how changing only the `SystemMessage` completely changes the model's tone and style, even though the question is identical.

## Statelessness

The model is stateless. This means that each call to `invoke()` is completely independent — the model has **no memory** of previous calls. Try this:

```python
response = llm.invoke("My name is John Doe")
print(response.text)

response = llm.invoke("What is my name?")
print(response.text)
```

The model won't remember your name, because the second call knows nothing about the first one. Each `invoke()` starts from scratch.

### Building a conversation manually

If we want the model to remember what was said before, **we** have to pass the full conversation history every time. This is where **messages** come in.

When you wrote `llm.invoke("What is the capital of France?")` earlier, LangChain silently wrapped your string into a `HumanMessage` — the same kind of object we're about to use explicitly. Under the hood, the model always works with a **list of messages**, not plain strings.

```python
from langchain_core.messages import HumanMessage

conversation_history = []

# First exchange
conversation_history.append(HumanMessage("My name is John Doe"))
response = llm.invoke(conversation_history)
conversation_history.append(response)

# Second exchange
conversation_history.append(HumanMessage("What is my name?"))
response = llm.invoke(conversation_history)
print(response.text)
```

Now the model sees the full conversation and can answer correctly. Notice that `response` (the model's reply) is itself a message object — an `AIMessage` — so we can append it directly to the list.

An immediate consequence of this behavior is that every time we call the model, **all** previous messages need to be passed again. This means the cost of each request grows as the conversation gets longer (more tokens in every request).

## Dissecting the structure of the response

Considering that we still have the response from the previous section, we can dissect the structure of the response.

```python
print("--- CONTENT ---")
print(resp.content)          # full structure returned by the API

# To extract just the text:
print(resp.text)

print("\n--- TOKENS ---")
print(resp.usage_metadata)        # {'input_tokens': ..., 'output_tokens': ..., 'total_tokens': ...}

print("\n--- STATUS/METADATA ---")
print(resp.response_metadata)     # e.g., {'finish_reason': 'STOP', ...}
```

> [!NOTE]
> With some models, `resp.content` returns a list of dictionaries like `[{'type': 'text', 'text': '...', 'extras': {'signature': '...'}}]` instead of a plain string. The `.text` property extracts just the text from the response (equivalent to accessing the `text` field inside the content structure). When you use `StrOutputParser()` in a chain, this extraction happens automatically.

As we can see, the output is not just the text of the response. It contains the content of the response, the usage metadata and the response metadata. Usage metadata is particularly important to understand the cost of the request.

### Understanding Tokens and Cost

LLM APIs charge based on **tokens** — roughly, pieces of words. A typical English word is 1–2 tokens. For example:
- "photosynthesis" ≈ 3 tokens
- "The cat sat on the mat" ≈ 7 tokens

Both input (your prompt) and output (the model's response) are counted. The usage metadata shows:
- `input_tokens`: tokens in your prompt
- `output_tokens`: tokens in the model's response
- `total_tokens`: sum of both

This matters because:
1. APIs charge per token (e.g., Gemini Flash: ~$0.075 per 1M input tokens)
2. Longer conversations cost more — you resend all history each time
3. Rate limits are often measured in tokens per minute

## The pipe operator in LangChain

LangChain provides a pipe operator that allows to chain multiple steps together. This is particularly useful to build a pipeline of steps. We saw that the response from a model is a complex AI message object. There could be cases in which we are interested in extracting only the content of the response. This can be easily achieved by composing dedicated functions with the pipe operator.

Let us consider an example in which we define a **chain** that summarizes a given text in a single sentence.

```python
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser

paragraph = (
    "CRISPR-Cas9 is a revolutionary genome-editing technology derived from a "
    "bacterial immune system. It uses a guide RNA to direct the Cas9 enzyme to "
    "a specific DNA sequence, where it creates a double-strand break. The cell's "
    "repair mechanisms then either disrupt the gene (via NHEJ) or insert new "
    "genetic material (via HDR)."
)

prompt = ChatPromptTemplate.from_template("You are a helpful assistant that summarizes text in a single sentence. {text}")

chain = prompt | llm | StrOutputParser()

result = chain.invoke({"text": paragraph})
print(result)
```

Alternatively, we can use the `from_messages` method to create a prompt template from a list of messages. This is useful in the case in which we want to use a system message to set the overall behavior of the model.

```python
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser

prompt = ChatPromptTemplate.from_messages([
    ("system", "You are a helpful assistant that summarizes text in a single sentence."),
    ("human", "{text}")
])

chain = prompt | llm | StrOutputParser()

result = chain.invoke({"text": paragraph})
print(result)
```

## Reusing prompts with `partial()` syntax

We can use the `partial` method to create a new prompt template from a base prompt template. This is useful to reuse the same prompt template with different parameters. In the example below we change the system prompt by adding a specific domain and a limit to the output.

```python
base = ChatPromptTemplate.from_messages([
    ("system", "You are a domain expert in {domain}. Keep output under {limit} words."),
    ("human", "{question}")
])

bio_prompt = base.partial(domain="bioinformatics", limit="60")
result = (bio_prompt | llm | StrOutputParser()).invoke({"question": "What is differential expression?"})
print(result)
```

The `partial()` method lets you **fill in some variables now** while leaving others for later. Think of it like a form with some fields already filled in — you can reuse the same base form with different pre-filled values (like `domain`) while leaving other fields (like `question`) for the user to complete each time.

## Conversation memory via `MessagesPlaceholder`

The `MessagesPlaceholder` is a special placeholder that allows to inject the conversation history into the prompt. This is useful to build chains where the context comes from the conversation history.

```python
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain_core.messages import HumanMessage, AIMessage

prompt = ChatPromptTemplate.from_messages([
    ("system", "You are an expert in pharmacology"),
    MessagesPlaceholder("history"),
    ("human", "{input}")
])

history = [
    HumanMessage(content="What is the putative target of venetoclax?"),
    AIMessage(content="The putative target of venetoclax is BCL2.")
]

chain = prompt | llm | StrOutputParser()
result = chain.invoke({"history": history, "input": "And what about Dabrafenib?"})
print(result)
```

> [!NOTE]
> The string `"history"` in `MessagesPlaceholder("history")` is a label you choose — it can be any name. When invoking the chain, you pass `{"history": [...]}` — the key must match the label you chose.

This template is particularly useful for injecting few-shot examples into the prompt.

## Swapping providers without refactoring the chains

As we have seen, defining a chain requires to specify the model provider and the model name. However, we can avoid this by defining a function that creates the chain for us. We just pass the provider and model name as arguments:

```python
# --- Configuration ---
# MODEL_NAME = "gemini-3-flash-preview"
# MODEL_PROVIDER = "google_genai"
MODEL_NAME = "my-model"
MODEL_PROVIDER = "openai"
BASE_URL = "https://myung-courageous-gino.ngrok-free.dev/v1"
API_KEY = "not-needed"

def make_chain(model_provider: str, model_name: str, temperature: float = 0.2, **kwargs):
    llm = init_chat_model(model=model_name,
                          model_provider=model_provider,
                          temperature=temperature,
                          **kwargs)
    prompt = ChatPromptTemplate.from_template("Give a one‑line TL;DR about {thing}.")
    return prompt | llm | StrOutputParser()

chain = make_chain(MODEL_PROVIDER, MODEL_NAME, base_url=BASE_URL, api_key=API_KEY)
result = chain.invoke({"thing": "G protein-coupled receptors"})
print(result)
```

## Exercises

The exercises for this lecture are in the `exercises/2_exercise/` folder. See [exercise_langchain.md](exercises/2_exercise/exercise_langchain.md) for instructions.

The exercises will guide you through:
1. **Exercise 1** (Easy): Build a chain that explains scientific concepts simply
2. **Exercise 2** (Easy): Simulate a doctor-patient conversation using message types
3. **Exercise 3** (Hard): Create a multi-specialist consultation system
