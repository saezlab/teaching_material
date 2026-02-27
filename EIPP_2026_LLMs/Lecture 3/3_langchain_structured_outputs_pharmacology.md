# LangChain Structured Outputs — Pharmacology Edition

> [!TIP]
> The code blocks in this file build on each other — each block assumes variables
> from the previous blocks exist. If you want to experiment, combine them into a
> single Python script or run them sequentially in a notebook.

## Table of Contents

1. [A pharmacology example](#1-a-pharmacology-example)
   - [Example 1: Drug monograph](#example-1-drug-monograph)
   - [Example 2: Oral vs injectable](#example-2-oral-vs-injectable)
2. [A nested example](#2-a-nested-example)
3. [A more general approach](#3-a-more-general-approach)

## Why pharmacology?

Pharmacology is a natural fit for practicing structured outputs because drug data is inherently structured. Every medication typically follows a predictable schema:
- **Identity**: Name and manufacturer.
- **Usage**: Indications and dosage.
- **Logistics**: Route of administration (oral, IV, etc.).

This mirrors the kind of structured information extraction you'll encounter in many biology research tasks, such as:
- Parsing **gene annotations**
- Extracting **experimental parameters** from papers
- Cataloguing **compound properties**

While the examples below focus on pharmacology, these patterns apply broadly across the life sciences.

---

> [!NOTE]
> The code examples in this file assume you have your environment set up from Lecture 2.
> The exercise scripts include `load_dotenv()` for you — in these reading examples, we skip
> it for brevity.

## 1. A pharmacology example

### Example 1: Drug monograph

Let us consider the case in which we want to get and parse the information from a drug monograph. Without structured outputs, we would have to manually parse the text and extract the information. This is error prone and requires a lot of code. With structured outputs, we can define a Pydantic model and have LangChain enforce the schema and validation for us. Let us first generate the drug monograph with a LLM:

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
drug = "aspirin"

generation_prompt = '''
<task>
Give me a paragraph with information about {drug}:  
 - name
 - indications
 - contraindications
 - common_adverse_effects
 - adult_dose_mg
 - routes
</task>
'''

response = llm.invoke(generation_prompt.format(drug=drug))
drug_monograph = response.content
```

Now that we have the response, let us parse it with a Pydantic model.

> [!NOTE]
> We could have used constrained generation to directly get the response in the correct format. However, here we break down the problem into smaller steps.
> Notably, this two-step process models the case in which we have knowledge from an external source in free text (say a publication) and we want to parse it into a structured format.

Let us now define a Pydantic model for the response.

> [!NOTE]
> **`Field()` patterns you'll see:**
> - `Field(...)` — required, no default (the `...` means "must be provided")
> - `Field(None, description="...")` — optional, defaults to `None`
> - `Field(description="...")` — required, with a description for documentation

```python
from pydantic import BaseModel, Field
from typing import List, Optional
from langchain.chat_models import init_chat_model

class Medication(BaseModel):
    name: str = Field(description="Generic or brand name")
    indications: List[str] = Field(description="Therapeutic indications")
    adult_dose_mg: Optional[float] = Field(None, description="Typical adult dose in mg if applicable")
    routes: List[str] = Field(description="Administration routes, e.g., oral, IV")
```

Now that we have the Pydantic model, we can use and ad hoc method provided by LangChain to constrained the model to generate a response abiding to the Pydantic model. Crucially, this method is called `with_structured_output` and it is available for all the models that support it (you can check [here](https://python.langchain.com/docs/integrations/chat/) the list of models that support it).

```python
structured_llm = llm.with_structured_output(Medication)
parsing_prompt = f"Extract medication info following the schema: {drug_monograph}"
monograph = structured_llm.invoke(parsing_prompt)
print(monograph)
```

> [!NOTE]
> Notice that `structured_llm.invoke()` returns a `Medication` object directly — not an
> `AIMessage` like the `.invoke()` calls in Lecture 2. The `with_structured_output()` method
> changes the return type to match your Pydantic model.

### Aside: Beyond pharmacology

The same pattern works for any domain. For example, to extract gene annotations:

```python
class GeneAnnotation(BaseModel):
    gene_symbol: str = Field(description="Official gene symbol, e.g., TP53")
    organism: str = Field(description="Species, e.g., Homo sapiens")
    function: str = Field(description="Brief functional description")

gene_llm = llm.with_structured_output(GeneAnnotation)
result = gene_llm.invoke("TP53 is a human tumor suppressor gene involved in cell cycle regulation and apoptosis.")
print(result)
```

The rest of this file focuses on pharmacology examples to build depth, but remember: these patterns apply to any structured data extraction task.

### Example 2: Oral vs injectable

Sometimes the same concept appears in different forms (oral vs injectable). A wrapper with optional fields lets the LLM choose the appropriate schema at runtime while keeping structure. However, since both fields are `Optional`, the LLM might fill in **both** — so we need two safeguards:
1. **Clear Field descriptions** that tell the LLM to only fill the matching field.
2. **A `@model_validator`** that enforces exactly one field is set, raising an error if the LLM gets it wrong.

Lets consider the following example:

> [!NOTE]
> **`@field_validator` vs `@model_validator`**
>
> - `@field_validator('field_name')` validates **one field** at a time. It uses `cls`
>   (the class) and `v` (the field value).
> - `@model_validator(mode="after")` validates the **entire model** after all fields
>   are set. It uses `self` (the full object), so you can check relationships between
>   fields — like ensuring exactly one of two optional fields is filled.

```python
from typing import Optional, List, Union
from pydantic import BaseModel, Field, model_validator

class OralMedication(BaseModel):
    name: str
    dosage_form: str = Field(description="e.g., tablet, capsule")
    dose_mg: Optional[float] = Field(None, description="Per dose in mg")

class InjectableMedication(BaseModel):
    name: str
    concentration_mg_per_ml: Optional[float] = Field(None, description="Strength in mg/mL")
    route: str = Field(description="e.g., IV, IM, SC")

class MedicationWrapper(BaseModel):
    """Wrapper that holds EITHER an oral OR an injectable medication, never both."""
    oral: Optional[OralMedication] = Field(None, description="Fill ONLY if the medication is oral")
    injectable: Optional[InjectableMedication] = Field(None, description="Fill ONLY if the medication is injectable")

    @model_validator(mode="after")
    def exactly_one_set(self):
        if self.oral is not None and self.injectable is not None:
            raise ValueError("Only one of oral or injectable should be set, not both")
        if self.oral is None and self.injectable is None:
            raise ValueError("At least one of oral or injectable must be set")
        return self

    def get_therapy_type(self):
        "Returns the type of the medication based on the non-None field"
        if self.oral is not None:
            return "oral"
        elif self.injectable is not None:
            return "injectable"
        else:
            return None


flexible_llm = llm.with_structured_output(MedicationWrapper)

oral_text = "Paracetamol 500 mg tablet for pain and fever (oral)."
inj_text = "Ondansetron injection 2 mg/mL, given IV for nausea."

prompt = "Extract the medication info. Fill ONLY the oral field if oral, or ONLY the injectable field if injectable: {text}"

oral = flexible_llm.invoke(prompt.format(text=oral_text))
inj = flexible_llm.invoke(prompt.format(text=inj_text))

print(oral.get_therapy_type())   # → oral
print(inj.get_therapy_type())    # → injectable
```


## 2. A nested example

Pharmacology often bundles related fields (dose + route + frequency). Nesting groups these logically. Add simple validators for safety. We keep it to three classes to stay readable.

```python
from pydantic import BaseModel, Field, field_validator
from typing import List, Optional

class DosageInstruction(BaseModel):
    route: str = Field(description="Route, e.g., oral, IV")
    dose_mg: float = Field(description="Dose per administration in mg")
    frequency: str = Field(description="e.g., q8h, once daily")

class Manufacturer(BaseModel):
    name: str = Field(description="Company name")
    country: Optional[str] = Field(None, description="Country of origin")

class Drug(BaseModel):
    name: str = Field(description="Generic name")
    indications: List[str] = Field(description="Therapeutic indications")
    dosage: DosageInstruction = Field(description="Standard adult dosing")
    manufacturer: Manufacturer = Field(description="Manufacturer details")

    @field_validator("dosage")
    def dose_positive(cls, v: DosageInstruction):
        if v.dose_mg <= 0:
            raise ValueError("dose_mg must be > 0")
        return v

# generation_prompt is defined in Example 1 above — now let's try with a different drug
drug = "ibuprofen"
drug_response = llm.invoke(generation_prompt.format(drug=drug))
drug_text = drug_response.content

nested_llm = llm.with_structured_output(Drug)
drug_profile = nested_llm.invoke(f"Extract detailed drug profile following the schema: {drug_text}")
```


## 3. A more general approach

The `with_structured_output` method is convenient, but not all LLM providers support it. The `PydanticOutputParser` offers an alternative that works with **any** model — it injects the schema directly into the prompt as text, so the model knows what structure to produce. This makes it a useful fallback when `with_structured_output` is unavailable.

```python
from langchain_core.output_parsers import PydanticOutputParser
parser = PydanticOutputParser(pydantic_object=Medication)
```

See how the parser gets initialized with the pydantic model we defined previously. Before wiring it into a chain, let's inspect what it actually tells the model to do:

```python
print(parser.get_format_instructions())
```

You'll see a block of text ending with a JSON schema inside backtick fences. To read it more clearly, pretty-print just the JSON part:

```python
import json

raw = parser.get_format_instructions()
json_str = raw.split("```")[1].strip()
schema = json.loads(json_str)
print(json.dumps(schema, indent=2))
```

Look at the output — you'll see the field names, types, and descriptions from your Pydantic model translated into a JSON schema. This is exactly what gets injected into the prompt to tell the model how to structure its response.

Now let's use these format instructions in a prompt template:

```python
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser

format_instructions = parser.get_format_instructions()

prompt = ChatPromptTemplate.from_messages([
    ("system", "Extract per schema:\n{format_instructions}"),
    ("human", "{text}"),
]).partial(format_instructions=format_instructions)

parsing_llm = prompt | llm | parser

# drug_text is the string we got from drug_response.content in Example 1
result = parsing_llm.invoke({"text": drug_text})
print(result)
```

## Exercises

The exercises for this section are in the `exercises/3_exercise/` folder. See [exercise_structured_outputs.md](exercises/3_exercise/exercise_structured_outputs.md) for instructions.

The exercises will guide you through:
1. **Exercise 1**: Use `Literal` to turn an LLM into a tissue classifier
2. **Exercise 2**: Define a Pydantic model for a dietary supplement and extract structured data
3. **Exercise 3**: Use optional fields with a model validator to classify route-specific therapies
4. **Exercise 4**: Build a nested model with validation for a vaccine profile