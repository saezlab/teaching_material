# Classical Programming vs. LLM Outputs: A Structural Perspective

The key insight of this lecture: LLMs are not just text generators — they can be **data producers**. By the end of this file, you'll understand why getting structured data out of LLMs matters and what challenges are involved.

## Table of Contents

1. [The Challenge of Unstructured LLM Outputs](#1-the-challenge-of-unstructured-llm-outputs)
2. [Structure in Classical Programming](#2-structure-in-classical-programming)
3. [Why Structure Matters in Production Systems](#3-why-structure-matters-in-production-systems)
4. [LLMs as Text Generators vs. Data Producers](#4-llms-as-text-generators-vs-data-producers)

## 1. The Challenge of Unstructured LLM Outputs

LLMs can perform for us a wide range of cognitive tasks that we might be interested in integrating in our applications. However, they naturally generate **free-form text** that humans understand but computers struggle to parse reliably. Consider the following example:

```python
from langchain.chat_models import init_chat_model

# --- Configuration ---
# MODEL_NAME = "gemini-3-flash-preview"
# MODEL_PROVIDER = "google_genai"
MODEL_NAME = "my-model"
MODEL_PROVIDER = "openai"
BASE_URL = "https://myung-courageous-gino.ngrok-free.dev/v1"
API_KEY = "not-needed"
llm = init_chat_model(model=MODEL_NAME, model_provider=MODEL_PROVIDER, base_url=BASE_URL, api_key=API_KEY, temperature=0.2)

# Ask for paper information
response = llm.invoke("""
Extract information about this paper:
"Attention Is All You Need" by Vaswani et al., published in 2017.
It introduced the Transformer architecture for neural networks.
""")

print(response.content)
```

After running the above code, a typical response from the LLM could look something like the following:

```text
This paper, titled "Attention Is All You Need," was authored by Vaswani and colleagues
in 2017. It's a groundbreaking work that introduced the Transformer architecture,
which revolutionized neural network design for natural language processing tasks.
The paper demonstrated that attention mechanisms alone, without recurrence or
convolution, could achieve state-of-the-art results in machine translation.
```

As we can see, the LLM has generated a text that is not structured in any way. This is a problem because it is not easy to parse this text into a structured format. Crucially we can pinpoint a few key issues:

* **Inconsistent format**: Each response might be structured differently
* **Missing fields**: No guarantee all required information is present
* **Extraction complexity**: Need to parse natural language to find data
* **Error-prone**: Typos, variations in formatting, incomplete information
* **No validation**: Can't programmatically verify data correctness


## 2. Structure in Classical Programming

In traditional software development, we work with **well-defined data structures** that provide predictability and type safety. For examples, let us consider the following class one could design to represent a research paper:

> [!NOTE]
> **Python Classes in 30 Seconds**
>
> A **class** is a blueprint for creating objects with specific fields. `__init__` is
> a special method that runs when you create a new object. `self` refers to the object
> being created. `self.title = title` stores the value so you can access it later as
> `paper.title`. If this is new to you, don't worry — Pydantic (File 2) will give you
> a much simpler way to define the same thing.

```python
# Classical approach: Structured data with clear types
class ResearchPaper:
    def __init__(self, title: str, authors: list[str], year: int, abstract: str):
        self.title = title
        self.authors = authors
        self.year = year
        self.abstract = abstract
        self.validate()

    def validate(self):
        if not self.title or len(self.title.strip()) == 0:
            raise ValueError("Title cannot be empty")
        if self.year < 1900 or self.year > 2026:
            raise ValueError("Year must be between 1900 and 2026")
        if not self.authors:
            raise ValueError("Must have at least one author")

# Usage is predictable and type-safe
paper = ResearchPaper(
    title="Attention Is All You Need",
    authors=["Vaswani, A.", "Shazeer, N.", "Parmar, N."],
    year=2017,
    abstract="The dominant sequence transduction models..."
)

# We can reliably access structured attributes
print(f"Paper from {paper.year}: {paper.title}")
print(f"Number of authors: {len(paper.authors)}")
```

As we can see this class has a number of **essential characteristics of classical structured data:**

* **Type safety**: Each field has a specific type
* **Validation**: Data constraints are enforced
* **Predictability**: You know exactly what fields exist
* **Error handling**: Clear failure modes when data is invalid


## 3. Why Structure Matters in Production Systems

As researchers, you will often want to use LLM outputs in your analysis scripts, databases, or pipelines. This only works if the data is **reliable and consistently structured**. Let us consider a few examples:

### Lab Database Integration

Imagine you are using an LLM to extract gene annotations from papers and store them in a database.

```python
# We want to store gene annotations extracted by the LLM
gene_db.insert({
    "gene_symbol": "???",    # How do we extract "TP53" reliably from free text?
    "organism": "???",        # What if the LLM says "humans" vs "Homo sapiens"?
    "function": "???",        # How do we split multiple functions?
    "chromosome": "???"       # What if it says "chromosome 17" vs just "17"?
})
```

### Automated Analysis Pipelines

Suppose you are screening drug compounds and need to compare properties extracted by the LLM.

```python
# Comparing IC50 values across compounds requires consistent structure
def find_most_potent(compound_list):
    for compound in compound_list:
        # These fail if the LLM returns inconsistent formats
        if compound.ic50_nm < 100:             # AttributeError if no .ic50_nm
            potent_hits.append(compound)

        if compound.target == "EGFR":          # What if it says "epidermal growth factor receptor"?
            egfr_compounds.append(compound)
```

### Batch Processing Research Papers

You might want to screen a large number of paper abstracts and filter them by topic or date.

```python
# Screening a batch of extracted paper metadata
def screen_recent_papers(paper_list):
    for paper in paper_list:
        # These fail with unstructured text
        if paper.year > 2020:                  # AttributeError if no .year
            recent_papers.append(paper)

        if "CRISPR" in paper.keywords:         # TypeError if keywords is a string
            crispr_papers.append(paper)
```

### Reproducible Experiment Logging

When logging experiment results extracted by an LLM, you need to know the data is complete and correctly formatted.

```python
# Logging experiment results needs consistent structure
def log_experiment(llm_response):
    try:
        result = parse_experiment_data(llm_response)
        validate_result(result)                # Were all required fields extracted?
        save_to_spreadsheet(result)            # Only works with structured data
    except ValidationError as e:
        print(f"Could not parse experiment data: {e}")
        # Without structure, we can't tell what went wrong
```

---

## 4. LLMs as Text Generators vs. Data Producers

### Traditional View: LLMs as Text Generators

The "traditional" use of LLMs (mostly through their web interface) is to use them as question-answering tools, in a similar fashion to how we use search engines:

```python
# LLM as a text completion tool
prompt = "Write a summary of the paper 'Attention Is All You Need'"
response = llm.invoke(prompt)
# Output: Natural language text for human consumption
```

### Modern View: LLMs as Data Producers/Integrators

However, the modern view of LLMs is that they are not only text generators, but also data producers and integrators. This is achieved by tasking them to produce outputs that can be parsed by computers in computer-readable formats.

```python
# LLM as a structured data extraction tool
prompt = """
Extract paper metadata in JSON format:
{
    "title": "exact title",
    "authors": ["author1", "author2"],
    "year": 2017,
    "abstract": "paper abstract"
}

Paper: "Attention Is All You Need" by Vaswani et al...
"""
response = llm.invoke(prompt)
# Expected output: Valid JSON that can be parsed programmatically
```

This shift in perspective opens up new possibilities and in the next section we will see how we can achieve this with bare prompt engineering.
