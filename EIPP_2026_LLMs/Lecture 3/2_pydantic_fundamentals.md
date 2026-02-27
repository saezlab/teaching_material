# Pydantic Fundamentals for Data Validation

> **A Note for Students**: Don't worry if you don't memorize every detail of this tutorial! Pydantic is a technical tool, and even experienced developers (and your instructors!) frequently check the documentation or use LLMs to get the formatting exactly right.
>
> The main goal today is to **understand what you can achieve** with Pydantic. Use this file as an entry point and a reference—something you can come back to whenever you need to build structured, reliable data for your LLM applications in the future.

> [!NOTE]
> **Installation:** Pydantic is installed automatically as a dependency of LangChain — no extra setup needed.
> If you want to make it explicit in your project: `uv add pydantic`

## Table of Contents

1. [Introduction to Pydantic](#1-introduction-to-pydantic)
2. [Basic Model Creation](#2-basic-model-creation)
   - [2.1 Essential Field Types](#21-essential-field-types)
   - [2.2 Creating Complex Models](#22-creating-complex-models)
3. [Custom Validation](#3-custom-validation)
4. [Field Constraints](#4-field-constraints)

## 1. Introduction to Pydantic

**Pydantic** is a Python library that uses **type annotations** to validate data and parse it into Python objects. It's particularly valuable for LLM applications because it provides a bridge between unstructured text and structured, validated data.

A Pydantic model is a Python class that inherits from BaseModel and defines a structured schema for data using type annotations. When you instantiate a Pydantic model, it automatically validates the input data against the declared types, coerces values where possible, and raises clear errors when the data doesn't conform. Essentially, it turns a plain class definition into a self-validating data container, giving you runtime type safety that standard Python type hints alone don't provide.

> **You can think of it as a contract**: you declare what shape your data should have (field names, types, constraints), and Pydantic enforces that contract every time data flows through the model. This makes it especially useful for parsing external inputs like API requests, configuration files, or database records, where you can't trust the data to be well formed.

### Examples

*   **A form you fill out**: Think of a Pydantic model like a structured form, say a registration form. The form has specific fields (name, email, age), each expecting a specific type of information. If you write your name where your age should go, the form rejects it. A Pydantic model works the same way: it defines the fields, the expected types, and refuses to accept data that doesn't fit.
*   **A recipe card**: A recipe card has a title (text), a prep time (a number in minutes), and a list of ingredients (a list of strings). If someone handed you a recipe card where the prep time said "banana," you'd know something is wrong. Pydantic does that check for you automatically: it looks at each field, verifies the data makes sense for that slot, and complains loudly if it doesn't.

### Why Pydantic for LLM Applications?

Consider the case in which we want to validate the data of a research paper we have extracted from a text. In order to do so we need to define a dedicated function that will validate the data for example:

```python
# Without Pydantic: Manual validation nightmare
def validate_paper_data(data):
    if not isinstance(data.get('title'), str):
        raise ValueError("Title must be a string")
    if not data['title'].strip():
        raise ValueError("Title cannot be empty")
    if not isinstance(data.get('year'), int):
        raise ValueError("Year must be an integer")
    if data['year'] < 1900 or data['year'] > 2026:
        raise ValueError("Year must be between 1900 and 2026")
    # ... many more validations

# Let's test this with some data
data = {"title": "Attention Is All You Need", "year": 2017}

# This works!
validate_paper_data(data)
print("Data is valid!")

# If we provided invalid data:
# invalid_data = {"title": "", "year": 1850}
# validate_paper_data(invalid_data)  <-- This would stop the program with a ValueError
```

By using Pydantic we can define a cleaner and more concise dedicated solution:

> [!NOTE]
> You may see two styles: `list[str]` (modern Python 3.9+) and `List[str]` (from the `typing` module).
> Both work identically with Pydantic. We use the `typing` import style throughout these exercises for consistency.

```python
# With Pydantic: Declarative and automatic
from pydantic import BaseModel, Field
from typing import List

class Paper(BaseModel):
    title: str = Field(..., min_length=1)         # string must be at least 1 character
    year: int = Field(..., ge=1900, le=2026)
    authors: List[str] = Field(..., min_length=1)  # list must contain at least 1 item

# Using the model with valid data
paper = Paper(
    title="Attention Is All You Need",
    year=2017,
    authors=["Ashish Vaswani", "Noam Shazeer"]
)

print(f"Paper created: {paper.title} ({paper.year})")

# Let's see what happens with invalid data:
try:
    bad_paper = Paper(title="", year=1850, authors=[])
except Exception as e:
    print(e)
    # Pydantic prints exactly which fields are invalid and why:
    # - title: String should have at least 1 character
    # - year: Input should be greater than or equal to 1900
    # - authors: List should have at least 1 item
```

> [!NOTE]
> The `...` (called "Ellipsis") is a special Python value. In Pydantic, `Field(...)` means
> "this field is **required** — you must provide a value." If you want a field to be optional,
> use `Field(default=None)` or `Optional[type] = None` instead.

This approach has several key benefits:
- Automatic validation: Type checking and constraint enforcement
- Clear error messages: Detailed feedback on validation failures
- JSON serialization: Easy conversion to/from JSON
- Documentation: Self-documenting schemas

**Why learn Pydantic in an LLM course?** In the next file, we'll use Pydantic models to tell LangChain *exactly* what structure we want from the LLM. The model definitions you learn here become the "contract" that the LLM must follow when generating structured data.

## 2. Basic Model Creation

### 2.1 Essential Field Types

When building Pydantic models, you'll work with various data types. Think of these types as the building blocks for your data structures. Let's explore each category with practical examples that you might encounter in real applications.

#### Basic Python Types

The foundation of any data model starts with simple, everyday data types:

```python
from pydantic import BaseModel

class ExperimentSample(BaseModel):
    # Text data — identifiers, labels, descriptions
    sample_id: str          # "CTRL-001", "DRUG-042"

    # Whole numbers — counts, replicates, passages
    replicate: int          # 1, 2, 3

    # Decimal numbers — measurements, concentrations
    concentration_uM: float # 10.5, 0.001, 100.0

    # True/False values — flags, conditions
    is_control: bool        # True, False

# Example: Creating an experiment sample
sample = ExperimentSample(
    sample_id="DRUG-042",
    replicate=3,
    concentration_uM=10.5,
    is_control=False
)
```

#### Collection Types for Multiple Values

When you need to store multiple related items, use collection types:

```python
from typing import List, Dict, Set, Tuple

class ResearchProject(BaseModel):
    # List: Ordered collection (can have duplicates)
    keywords: List[str]              # ["AI", "machine learning", "neural networks"]

    # Dictionary: Key-value pairs for structured data
    metadata: Dict[str, str]         # {"funding": "NSF", "status": "active"}

    # Set: Unique values only (no duplicates)
    collaborators: Set[str]          # {"Alice", "Bob", "Charlie"}

    # Tuple: Fixed-size, ordered collection
    coordinates: Tuple[float, float] # (40.7128, -74.0060) for lat/lng

# Example: Defining a research project
project = ResearchProject(
    keywords=["deep learning", "computer vision", "medical imaging"],
    metadata={"grant_id": "NSF-2024-001", "duration": "3 years"},
    collaborators={"Dr. Smith", "Dr. Johnson", "Dr. Chen"},
    coordinates=(42.3601, -71.0589)  # Boston coordinates
)
```

> [!NOTE]
> `Set` and `Tuple` are shown here for completeness. The exercises in this lecture focus
> on `List`, `Dict`, and `Optional` — you won't need `Set` or `Tuple` for the exercises.

#### Flexible Types for Variable Data

Sometimes your data can be different types or might not always be present:

```python
from typing import Optional, Union

class FlexibleUserData(BaseModel):
    # Optional: Field might be missing or None
    middle_name: Optional[str] = None    # Can be "Marie" or None

    # Union: Field can be one of several types
    user_id: Union[str, int]            # Can be "user_123" or 12345

    # Combined: Optional field that accepts multiple types
    phone: Optional[Union[str, int]] = None  # "555-1234", 5551234, or None

# Examples showing flexibility
user1 = FlexibleUserData(user_id="USER001")  # String ID, no middle name or phone
user2 = FlexibleUserData(
    user_id=42,
    middle_name="Elizabeth",
    phone="555-0123"
)
```

Pydantic supports many more types (dates, UUIDs, enums, etc.) — see the [Pydantic documentation](https://docs.pydantic.dev/) for details. For this course, the types above are all you need.

#### Constrained Values with Literal

Sometimes a field should only accept specific, predefined values — for example, a fixed set of categories. Use `Literal` for this:

```python
from typing import Literal
from pydantic import BaseModel

class ExperimentStatus(BaseModel):
    # Literal: Field MUST be one of the specified values — anything else is rejected
    status: Literal["running", "completed", "failed"]
    organism: Literal["human", "mouse", "rat", "zebrafish"]

# Valid:
exp = ExperimentStatus(status="running", organism="mouse")
print(exp)

# Invalid — Pydantic will raise a ValidationError:
# exp = ExperimentStatus(status="paused", organism="dog")
```

`Literal` is especially useful when you want an LLM to **classify** something into a fixed set of categories — you'll use this pattern in Exercise 1.

### 2.2 Creating Complex Models

Pydantic allows also to create more complex models by nesting models inside other models. Let us see an example using again the case of a research paper. Let us start by defining a simple model for an author:

```python
from pydantic import BaseModel
from typing import List, Optional
from datetime import datetime

class Author(BaseModel):
    name: str
    email: Optional[str] = None
    affiliation: Optional[str] = None
```

We can now create a more complex model for a research paper by nesting the `Author` model inside the `ResearchPaper` model:

```python
class ResearchPaper(BaseModel):
    title: str
    authors: List[Author]
    year: int
    abstract: str
    keywords: List[str]
    published_date: Optional[datetime] = None
```

Let us now instantiate both models:

```python
# Creating instances
author1 = Author(name="John Doe", email="john@university.edu")
author2 = Author(name="Jane Smith", affiliation="MIT")

paper = ResearchPaper(
    title="Attention Is All You Need",
    authors=[author1, author2],
    year=2017,
    abstract="This paper introduces the Transformer...",
    keywords=["transformer", "attention", "neural networks"]
)

print(paper.title)  # Access fields normally
print(paper.model_dump())  # Convert to dictionary
print(paper.model_dump_json())  # Convert to JSON string
```

> [!NOTE]
> **What are decorators and validators?**
>
> A **decorator** is a line starting with `@` placed above a function. It modifies how
> that function behaves. In Pydantic, `@field_validator('title')` means "run this
> function automatically every time someone sets the `title` field."
>
> The validator function receives two arguments:
> - `cls` — the class itself (required by Pydantic, similar to `self`)
> - `v` — the value that was passed for that field
>
> The validator must return a value: either `v` unchanged (if valid), a modified version
> (to clean or format the data), or raise a `ValueError` to reject the input.

## 3. Custom Validation

Pydantic allows also to define custom validation rules for the fields of the model. Let us see an example:

```python
from pydantic import BaseModel, field_validator

class Paper(BaseModel):
    title: str
    year: int
    journal: str
    impact_factor: float

    @field_validator('title')
    def title_must_not_be_empty(cls, v):
        if not v.strip():
            raise ValueError('Title cannot be empty')
        return v.strip().title()  # Clean and format

    @field_validator('year')
    def year_must_be_reasonable(cls, v):
        if v < 1900 or v > 2026:
            raise ValueError('Year must be between 1900 and 2026')
        return v

    @field_validator('impact_factor')
    def impact_factor_positive(cls, v):
        if v < 0:
            raise ValueError('Impact factor must be positive')
        return v

    @field_validator('journal')
    def journal_format(cls, v):
        # Standardize journal names
        journal_mappings = {
            'nature': 'Nature',
            'science': 'Science',
            'cell': 'Cell'
        }
        return journal_mappings.get(v.lower(), v)

```

Validators are defined in the class statement and are defined as methods with the `@field_validator` decorator. The method receives the value of the field as input and returns the validated value.

```python
# Usage
paper = Paper(
    title="  attention is all you need  ",  # Will be cleaned and formatted
    year=2017,
    journal="nature",  # Will be standardized to "Nature"
    impact_factor=42.8
)

print(paper.title)    # "Attention Is All You Need"
print(paper.journal)  # "Nature"
```

## 4. Field Constraints

Certain types of constraints can also be defined using the `Field` function. Let us see an example:

```python
from pydantic import BaseModel, Field
from typing import List

class ConstrainedPaper(BaseModel):
    # String constraints
    title: str = Field(..., min_length=5, max_length=200)
    abstract: str = Field(..., min_length=50)

    # Numeric constraints
    year: int = Field(..., ge=1900, le=2026)  # ge = greater or equal, le = less or equal
    pages: int = Field(..., gt=0, lt=1000)    # gt = greater than, lt = less than
    impact_factor: float = Field(..., ge=0.0)

    # Collection constraints
    authors: List[str] = Field(..., min_length=1, max_length=20)
    keywords: List[str] = Field(default=[], max_length=10)
```

For more advanced use cases, `Field` also supports descriptions, examples, and regex pattern validation:

```python
from pydantic import BaseModel, Field

class DOIRecord(BaseModel):
    title: str = Field(..., min_length=5)
    doi: str = Field(
        ...,
        description="Digital Object Identifier",
        examples=["10.1000/182"],
        pattern=r"^10\.\d{4,}/.*"  # DOI format validation
    )
```

```python
# If we try to create a paper that breaks these rules:
# paper = ConstrainedPaper(
#     title="AI",               # Error: Too short
#     abstract="Short abstract", # Error: Too short
#     year=1800,               # Error: Too old
#     pages=0,                 # Error: Must be > 0
#     authors=[],              # Error: Empty list
# )
# Pydantic will stop the program and show us exactly which rules were broken!
```

## Coming up next

In the next part of the lecture we will see how to use Pydantic to validate data from LLMs.
