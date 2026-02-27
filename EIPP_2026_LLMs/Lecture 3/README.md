# Lecture 3: Structured Outputs from LLMs

## Overview

This lecture focuses on extracting structured, reliable data from Large Language Models (LLMs). While LLMs naturally generate unstructured text, modern applications often require structured outputs like JSON objects, validated data schemas, or specific formats. Students will learn how to move from free-form text generation to predictable, parseable outputs that can be integrated into larger software systems.

## Learning Objectives

By the end of this lecture, students will be able to:

- Understand the importance of structured outputs in production LLM applications
- Compare classical structured programming approaches with LLM text generation
- Extract JSON data from LLM responses using manual parsing techniques
- Create and validate data schemas using Pydantic models
- Use LangChain's structured output features for reliable data extraction
- Build robust pipelines that combine LLM intelligence with structured data requirements

## Lecture Materials Index

### 1. Classical Programming vs. LLM Outputs
**📁 File:** [`0_classical_vs_llm_structure.md`](0_classical_vs_llm_structure.md)

Reflect on how traditional software handles structured data and compare it with the challenges of extracting structure from LLM-generated text.

**Key Topics:**
- Structured data in classical programming
- Challenges with unstructured LLM outputs
- Why structure matters in production systems
- Trade-offs between flexibility and reliability

---

### 2. Manual JSON Parsing from LLMs
**📁 File:** [`1_manual_prompting.md`](1_manual_prompting.md)

Learn foundational techniques for extracting JSON data from LLM responses, including prompt engineering and robust parsing strategies.

**Key Topics:**
- Prompting for JSON output
- Handling malformed JSON responses
- Error handling and retry strategies
- Performance considerations

---

### 3. Introduction to Pydantic Models
**📁 File:** [`2_pydantic_fundamentals.md`](2_pydantic_fundamentals.md)

Master Pydantic for data validation and schema definition, creating the foundation for reliable structured outputs.

**Key Topics:**
- Pydantic model basics and syntax
- Type annotations and validation
- Setting rules and limits on fields
- Nested models and complex structures
- Error handling and validation messages

---

### 4. LangChain Structured Outputs — Pharmacology Edition
**📁 File:** [`3_langchain_structured_outputs_pharmacology.md`](3_langchain_structured_outputs_pharmacology.md)

Leverage LangChain's structured output features with domain examples from pharmacology. You'll combine LLM generation with Pydantic validation and see unions and nested schemas in action.

**Key Topics:**
- `with_structured_output()` method
- Pydantic integration in LangChain
- Handling different forms of data (e.g., oral vs injectable medications)
- Building models within models (nesting)
- An alternative parsing approach for broader model compatibility

---

### 5. Exercises
**📁 Folder:** [`exercises/3_exercise/`](exercises/3_exercise/exercise_structured_outputs.md)

The exercises span concepts from all four files above — from Pydantic basics to LangChain integration.

---

## Prerequisites

- Completion of Lecture 2 (LLM APIs and LangChain)
- Python programming experience
- Basic understanding of JSON and data structures
- Familiarity with type hints (helpful but not required)


## Getting Started

> [!NOTE]
> No new package setup is needed before reading these files. When you reach the exercises,
> the exercise instructions will guide you through any additional setup.

1. **Foundation:** Begin with [`0_classical_vs_llm_structure.md`](0_classical_vs_llm_structure.md) to understand the conceptual framework
2. **Manual Techniques:** Work through [`1_manual_prompting.md`](1_manual_prompting.md) for foundational skills
3. **Schema Design:** Master [`2_pydantic_fundamentals.md`](2_pydantic_fundamentals.md) for robust data validation
4. **Framework Integration:** Apply [`3_langchain_structured_outputs_pharmacology.md`](3_langchain_structured_outputs_pharmacology.md) for production-ready solutions
5. **Exercises:** Complete the exercises in [`exercises/3_exercise/`](exercises/3_exercise/exercise_structured_outputs.md)

## Practical Applications

This lecture prepares you for real-world scenarios such as:

- **Research Data Extraction:** Automatically extracting structured metadata from scientific papers
- **Medical Information Processing:** Converting clinical notes into structured medical records
- **API Integration:** Creating reliable interfaces between LLMs and bioinformatics pipelines
- **Clinical Trial Data Processing:** Structuring trial protocols, endpoints, and patient data from unstructured sources
- **Literature Mining:** Extracting gene-disease associations, drug targets, or experimental conditions from publications

## Note on Model Capabilities

Different LLM providers have varying levels of support for structured outputs. Some models (like newer OpenAI models) have native structured output modes, while others require careful prompt engineering. This lecture covers techniques that work across different providers.

## Course Context

This lecture is part of **Topic 1: Data Integration and LLMs** in the EIPP 2026 LLMs. It bridges the gap between basic LLM usage (Lecture 2) and advanced data integration techniques that will be covered in subsequent lectures.

## Additional Resources

- [Pydantic Documentation](https://docs.pydantic.dev/)
- [LangChain Structured Outputs](https://python.langchain.com/docs/how_to/structured_output/)
- [JSON Schema Specification](https://json-schema.org/)
- [OpenAI Structured Outputs Guide](https://platform.openai.com/docs/guides/structured-outputs)
