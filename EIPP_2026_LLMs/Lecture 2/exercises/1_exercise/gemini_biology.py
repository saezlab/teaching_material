"""
Gemini Biology Exercise
========================
Sends a biological text to the Gemini API and asks it to
extract and explain the key concepts.

Usage (from the playground/ root):
    uv run lecture_2/1_gemini_api/gemini_biology.py
"""

import os
import sys

try:
    from google import genai
except ImportError:
    print("Error: the 'google-genai' package is not installed.", flush=True)
    print("Run this from the playground/ root:  uv add google-genai", flush=True)
    sys.exit(1)

BIOLOGICAL_TEXT = """\
The CRISPR-Cas9 system, originally discovered as an adaptive immune mechanism \
in prokaryotes, has been repurposed as a powerful genome-editing tool. The system \
relies on a guide RNA (gRNA) that directs the Cas9 endonuclease to a specific \
genomic locus through complementary base pairing. Upon binding, Cas9 introduces \
a double-strand break (DSB) in the target DNA. The cell then repairs the break \
through one of two major pathways: non-homologous end joining (NHEJ), which \
often introduces insertions or deletions (indels), or homology-directed repair \
(HDR), which can incorporate a donor template for precise editing. Recent \
advances include base editors and prime editors, which enable targeted nucleotide \
changes without requiring double-strand breaks.\
"""

PROMPT = f"""\
You are a molecular biology teaching assistant.

A student gives you the following text:

---
{BIOLOGICAL_TEXT}
---

Please:
1. List the key biological concepts mentioned in the text.
2. Briefly explain each concept in one or two sentences, \
suitable for an undergraduate student.
3. Suggest one follow-up question the student could explore.
"""


def main():
    api_key = os.environ.get("GOOGLE_API_KEY")
    if not api_key:
        print("Error: GOOGLE_API_KEY environment variable is not set.", flush=True)
        print("Run: export GOOGLE_API_KEY='your-api-key-here'", flush=True)
        sys.exit(1)

    client = genai.Client(api_key=api_key)

    print("Sending biological text to Gemini...\n", flush=True)

    try:
        response = client.models.generate_content(
            model="gemini-3-flash-preview",
            contents=PROMPT,
        )
    except Exception as e:
        error_text = str(e)
        if "429" in error_text or "RESOURCE_EXHAUSTED" in error_text:
            print("Rate limit hit (429). The free Gemini API has strict quotas.")
            print('Wait about 60 seconds and try again, or change the model to "gemini-2.5-flash" in this script.')
            sys.exit(1)
        raise

    print(response.text)


if __name__ == "__main__":
    main()
