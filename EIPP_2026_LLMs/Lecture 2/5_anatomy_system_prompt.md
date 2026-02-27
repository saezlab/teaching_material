# (Optional) Anatomy of a System Prompt

## Why This Section Exists

So far, we've used simple system prompts like "You are a helpful assistant" or "You are a biology tutor." But in real applications — especially agents that use tools and need to follow complex workflows — system prompts become much more detailed and structured.

The **Biomni agent** (a biology research assistant built at Stanford) is a good example of a **production-grade system prompt**. By studying it, you'll see how professional LLM applications structure their instructions to ensure consistent, reliable behavior.

**This section is optional**, but highly recommended if you plan to build your own agents in later lectures.

## Key Components of a System Prompt

A well-designed system prompt typically includes:

1. **Role Definition** — Who is the agent? What is its purpose?
   - Example: "You are Biomni, a biology research assistant specializing in molecular biology and drug discovery."

2. **Behavioral Guidelines** — How should the agent respond?
   - Example: "Always cite sources," "Explain technical concepts in simple terms," "Admit when you don't know something."

3. **Available Tools** (for agents) — What actions can the agent take?
   - Example: "You have access to 50+ tools including PubMed search, protein structure lookup, and pathway analysis."

4. **Output Format** — How should responses be structured?
   - Example: "Format tool calls as JSON," "Use markdown for formatting," "Limit responses to 200 words."

5. **Constraints and Rules** — What should the agent NOT do?
   - Example: "Never make up data," "Don't perform medical diagnosis," "Stay within the biology domain."

## The Biomni Prompt

The Biomni agent uses a ~2-page system prompt that exemplifies all the components above. You can find it in the [Biomni paper](https://biomni.stanford.edu/paper.pdf) (search for "Biomni prompt" — it appears in the supplementary materials).

**What to look for when reading it:**
- Notice how the role is defined in the first paragraph
- Count how many tools are listed and how they're described
- Identify the output format requirements (JSON structure, citation format)
- Look for constraints — statements with "always," "never," or "must"

## Comparing to Our Simple Prompts

Here's how the Biomni prompt compares to what we've used so far:

| Our prompts | Biomni prompt |
|-------------|---------------|
| 1 sentence | ~2 pages |
| "You are a helpful biology tutor" | Detailed role with specialty areas |
| No tools | 50+ tools with descriptions |
| No format requirements | Strict JSON output format |
| No explicit constraints | 10+ behavioral rules |

**When do you need a complex prompt like Biomni's?**
- When building agents that use tools (Lecture 4–5)
- When consistency across many requests is critical
- When working in regulated domains (medical, legal)
- When you need to constrain model behavior tightly

**When is a simple prompt enough?**
- One-off queries or simple Q&A
- Exploratory data analysis
- Non-critical applications where variability is acceptable

## Try It Yourself (Optional Exercise)

Draft a system prompt for a **"Protein Structure Tutor"** agent that includes:
1. **Role:** A tutor for graduate students learning protein structural biology
2. **Behavior:** Patient, uses analogies, cites PDB entries when relevant
3. **Tools (hypothetical):** PDB structure lookup, 3D visualization, sequence alignment
4. **Output format:** Markdown with headers, bullet points for key facts
5. **Constraints:** Never makes up PDB IDs, admits uncertainty, stays within structural biology

This exercise will prepare you for Lecture 4–5 where we build agents with real tools.
