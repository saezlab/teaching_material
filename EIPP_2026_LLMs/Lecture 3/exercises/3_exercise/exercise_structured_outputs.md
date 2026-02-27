# Exercise — Structured Outputs from LLMs

## Goal

Practice using Pydantic models with LangChain's structured output features to extract validated, structured data from LLM responses.

## Prerequisites

- You have completed the Lecture 2 exercises (playground setup, Gemini API, LangChain basics)

## Setup

### 1. Install Pydantic (if not already installed)

Pydantic is already installed as a dependency of LangChain, but let's make it explicit in our project:

> [!TIP] Terminal
> ```bash
> uv add pydantic
> ```

> [!NOTE]
> The exercise scripts use the EBI GPU cluster endpoint (configured in Lecture 2). If the
> endpoint is unavailable, uncomment the Gemini lines in the configuration block at the
> top of each script.

> [!NOTE]
> **Pydantic warnings you can ignore:** When running exercises that use `with_structured_output()`,
> you may see warnings like `PydanticSerializationUnexpectedValue`. These are harmless —
> they come from an internal interaction between LangChain and Pydantic. Your code is working
> correctly. The starter scripts already suppress these warnings for you.

> [!NOTE]
> **Speed:** Exercises that make two LLM calls (generate text + parse it) may take a minute
> or two to complete. This is normal — just wait for the output.

### 2. Create the exercise folder and copy the starter scripts

> [!TIP] Terminal
> ```bash
> mkdir -p lecture_3/3_structured_outputs
> cp "../Lecture 3/exercises/3_exercise/verify_setup.py" lecture_3/3_structured_outputs/
> cp "../Lecture 3/exercises/3_exercise/exercise_1_tissue_classifier.py" lecture_3/3_structured_outputs/
> cp "../Lecture 3/exercises/3_exercise/exercise_2_supplement.py" lecture_3/3_structured_outputs/
> cp "../Lecture 3/exercises/3_exercise/exercise_3_union_classification.py" lecture_3/3_structured_outputs/
> cp "../Lecture 3/exercises/3_exercise/exercise_4_nested_vaccine.py" lecture_3/3_structured_outputs/
> ```

> [!NOTE]
> You can also copy the files using your operating system's file explorer instead of the terminal. Just drag or copy the four Python files from `Lecture 3/exercises/3_exercise/` into `playground/lecture_3/3_structured_outputs/`.

All commands below assume you are in the `playground/` folder.

### 3. Verify your setup

> [!TIP] Terminal
> ```bash
> uv run lecture_3/3_structured_outputs/verify_setup.py
> ```

You should see `[OK]` for each check and "All checks passed" at the end. If anything fails, follow the instructions in the error message before continuing.

---

## Exercise 1: LLM as a Tissue Classifier (Easy)

In this exercise you will use a Pydantic model with a `Literal` field to turn an LLM into a **classifier**. Given five histological text excerpts, the LLM must pick the organ of origin from a fixed set of options.

This is a powerful pattern: by constraining the output to a `Literal`, the LLM can only return one of the allowed values — turning free-text generation into structured classification.

**Your tasks** (edit `lecture_3/3_structured_outputs/exercise_1_tissue_classifier.py`):

1. Define a `TissueClassification` Pydantic model with:
   - `organ: Literal["heart", "liver", "brain", "bone", "skin"]` — the classified organ
   - `key_features: List[str]` — histological features found in the excerpt
2. Create a structured LLM with `llm.with_structured_output(TissueClassification)`
3. Loop over the five provided excerpts, classify each one, and print the results

> [!TIP] Terminal
> ```bash
> uv run lecture_3/3_structured_outputs/exercise_1_tissue_classifier.py
> ```

**Expected outcome:** Each excerpt is correctly classified (heart, liver, brain, bone, skin) and the key features match the histological terms in the text.

---

## Exercise 2: Structured Output for a Dietary Supplement (Easy)

In this exercise you will define a Pydantic model for a dietary supplement and use LangChain's `with_structured_output()` to extract structured data from LLM-generated text.

**Your tasks** (edit `lecture_3/3_structured_outputs/exercise_2_supplement.py`):

1. Define a `Supplement` Pydantic model with these fields:
   - `common_names: List[str]` — common names for the supplement
   - `active_components: List[str]` — active chemical components
   - `uses: List[str]` — therapeutic or health uses
   - `typical_adult_dose_mg: Optional[float]` — typical adult dose in milligrams
   - `forms: List[str]` — available forms (e.g., tablet, capsule, liquid)
2. Use the LLM to generate a short paragraph about "vitamin D3" covering all these fields
3. Use `with_structured_output(Supplement)` to parse the generated text into a `Supplement` object
4. Print and verify the parsed object (check that the types and values look reasonable)

> [!TIP] Terminal
> ```bash
> uv run lecture_3/3_structured_outputs/exercise_2_supplement.py
> ```

**Expected outcome:** A `Supplement` object with populated fields matching the schema. You should see the common names, active components, uses, dose, and forms printed out.

---

## Exercise 3: Union Classification for Route-Specific Therapies (Medium)

In this exercise you will define two Pydantic models for different therapy routes and use a **wrapper with optional fields** so the LLM correctly classifies the therapy type — the same pattern from Example 2 in the lecture.

**Your tasks** (edit `lecture_3/3_structured_outputs/exercise_3_union_classification.py`):

1. Define `TopicalTherapy` with these fields:
   - `name: str` — name of the topical therapy
   - `formulation: str` — e.g., cream, ointment, gel
   - `site_of_application: str` — where it is applied
2. Define `InhaledTherapy` with these fields:
   - `name: str` — name of the inhaled therapy
   - `device_type: str` — e.g., HFA inhaler, DPI
   - `dose_mcg_per_puff: Optional[float]` — dose per puff in micrograms
3. Create a `TherapyWrapper` model with two optional fields:
   - `topical: Optional[TopicalTherapy]` — fill only if the therapy is topical
   - `inhaled: Optional[InhaledTherapy]` — fill only if the therapy is inhaled
4. Add a `@model_validator` that enforces exactly one of the two fields is set (not both, not neither)
5. Parse two inputs: a topical therapy text and an inhaled therapy text
6. Print the type detected for each

**Why a model validator?** Since both fields are `Optional`, the LLM might fill in both — e.g. populating `TopicalTherapy` fields for an inhaled therapy. The `@model_validator` catches this and raises an error, ensuring the wrapper always holds exactly one therapy type.

> [!TIP] Terminal
> ```bash
> uv run lecture_3/3_structured_outputs/exercise_3_union_classification.py
> ```

**Expected outcome:** The script correctly classifies "topical" for the cream input and "inhaled" for the inhaler input. The validator ensures exactly one field is populated.

---

## Exercise 4: Nested Model with Validation — Vaccine Profile (Hard)

In this exercise you will build nested Pydantic models with a field validator for a vaccine profile.

**Your tasks** (edit `lecture_3/3_structured_outputs/exercise_4_nested_vaccine.py`):

1. Define `DosePlan` with these fields:
   - `route: str` — administration route (e.g., intramuscular)
   - `dose_value: float` — the dose amount
   - `dose_unit: str` — the unit (e.g., mL, mg)
2. Define `ManufacturerInfo` with these fields:
   - `name: str` — manufacturer company name
   - `country: Optional[str]` — country of origin
3. Define `VaccineProfile` with these fields:
   - `name: str` — vaccine name
   - `indications: List[str]` — what the vaccine protects against
   - `schedule: DosePlan` — dosing information (nested model)
   - `manufacturer: ManufacturerInfo` — manufacturer details (nested model)
4. Add a `@field_validator` on `VaccineProfile` ensuring `dose_value > 0`
5. Use the LLM to generate a paragraph about a vaccine, then parse it with `with_structured_output(VaccineProfile)`
6. Print the validated object and its nested fields

> [!TIP] Terminal
> ```bash
> uv run lecture_3/3_structured_outputs/exercise_4_nested_vaccine.py
> ```

**Expected outcome:** A `VaccineProfile` object with nested `DosePlan` and `ManufacturerInfo` objects, all properly validated.

---

## What your playground should look like now

```
playground/
├── .env
├── .gitignore
├── .venv/
├── pyproject.toml
├── uv.lock
├── lecture_2/
│   ├── 0_uv_setup/
│   ├── 1_gemini_api/
│   ├── 2_langchain/
│   ├── 3_openrouter/
│   └── 4_ollama/
└── lecture_3/
    └── 3_structured_outputs/
        ├── verify_setup.py
        ├── exercise_1_tissue_classifier.py
        ├── exercise_2_supplement.py
        ├── exercise_3_union_classification.py
        └── exercise_4_nested_vaccine.py
```
