# Topic 1: Data Integration and LLMs

## Trainer

Francesco Carli, EMBL-EBI & OpenTargets, UK

## Project overview

In this project, we will explore how LLMs can be applied to real-world
biological problems through a hands-on, hackathon-style approach. Starting from
a concrete biological question, we will guide participants in building a
functional prototype that leverages LLMs to interpret and connect diverse data
modalities. Key topics will include strategies to enhance model factuality and
knowledge access, such as vector-based retrieval and prompt-to-query methods. We
will also demonstrate how LLMs can interact with external tools to improve their
problem-solving capabilities and introduce techniques for generating structured
outputs. Crucially, evaluation methods for LLMs will be covered, with a focus on
relevance, coherence, and scientific utility. The project emphasizes a
practical, iterative mindset, encouraging participants to test, refine, and
critically assess their solutions. By the end, attendees will gain both
conceptual insight and hands-on experience in using LLMs as adaptable components
within biological data analysis workflows.

## Project structure

We will start by a short recap of the most important technologies (by Francesco), [Lecture 1](#lecture-1-theory). After this, we will organise ourselves into groups and come up with pragmatic aims for the week for each group. These can be informed by your research interests or general likings. Importantly, they should have the right scope to be tackled alongside the other lessons during the rest of the week, to be presented on Friday.

The other lessons will be dynamically integrated into the workflow of approaching each group's goals. We will switch between general explanations and problem-driven responsive lectures. The topics we aim to cover during the week are:

- [Inference](#lecture-2-inference)
- [Structured Outputs](#lecture-3-structured-outputs)
- [Tool Usage](#lecture-4-tool-usage)
- [Agents](#lecture-5-agents)

## The playground

The `playground/` folder is the **shared Python environment** for all exercises across all lectures. Students set it up once in Exercise 0 of Lecture 2 and keep building on it throughout the course. It is managed by **UV** with Python 3.11. Exercise work is organized inside it by lecture (e.g., `playground/lecture_2/0_uv_setup/`).

### Lecture 1 ("theory")

[Content](Lecture%201)

- Next token prediction
- Pre-training
- Instruction tuning (RLHF)
- Post-training (math and tool usage)
- Mixtures-of-experts (MOE)
- Thinking models vs standard models
- Commercial vs "open source models"
- Where to find open models
- Model dimension
- Quantisation
- Where to find closed models
- Some questions/suggestions to reflect on possible projects

### Lecture 2 (Inference)

[Content](Lecture%202)

- Setting up the environment with UV and creating the playground
- Getting a free API key for Gemini
- Introduction to LangChain and how it can be used to call models
	* Setting environment variables for API keys
	* `init_chat_model` for provider-agnostic model creation
	* Human, AI and system messages (conversation templates)
	* The `invoke` method and dissecting its output (content, tokens, status)
	* Chains and the pipe operator
	* Prompt templates and `MessagesPlaceholder` for conversation memory
- Extending LangChain to use OpenRouter models
	* Getting an OpenRouter API key
	* Using `ChatOpenAI` with OpenRouter's OpenAI-compatible API
	* Calling free models from OpenRouter
- Local inference with Ollama
	* Installing and running Ollama
	* Pulling and running models locally (`qwen3:0.6b`)
	* Integrating Ollama with LangChain
	* Comparing local vs cloud inference
- (Optional) Anatomy of a system prompt (the Biomni example)

### Lecture 3 (Structured outputs)

[Content](Lecture%203)

- Reflection on how classical code is structured
- Have the model write an output in JSON format and write a dedicated parser to extract it
- Introduction to Pydantic models (a few examples, show types, and the field decorator)
- Show how to use LangChain to get structured outputs automatically:
	* `with_structured_output()` method
	* `PydanticOutputParser` (more general approach)

### Lecture 4 (Tool usage)

[Content](Lecture%204)

- Some reference on tools (https://arxiv.org/pdf/2305.15334)
- Simple example with the weather (the classical one)
- Exposing a tool to the LLM from a conceptual point of view (docstring and signature)
- Trying to call a function from scratch (explain to the model what the tool does and how it should call it. Parse the output and try to call the tool. We have an example prompt in BioChatter)
- A tool in LangChain + tool binding
- Invoking the tool with the result from the model
- An example of a helpful tool (re-use the example from BioChatter);
- Introduction to MCP
	* What is its purpose
	* how can be used in Python
- How to attach an MCP server to LangChain (mcp adapter libraries)
- A simple example of a local MCP server (calculator)
- A simple example of a remote MCP (biocontext.AI)

### Lecture 5 (Agents)

[Content](Lecture%205)

- A few words on agents
- Brief introduction to the ReAct framework
- Brief introduction to langgraph
	* how to define nodes and edges
	* how to define a conditional edge
	* tool node
	* graph state
	* invoke the graph
- Replicate ReAct in langgraph
