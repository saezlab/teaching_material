"""
Solution — Exercise 3: Multi-Specialist Consultation System

This script makes 3 API calls total:
  - 2 in step 4 (one per specialist)
  - 1 in step 5 (follow-up with cardiology)
"""

from dotenv import load_dotenv
load_dotenv()

from langchain.chat_models import init_chat_model
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain_core.messages import HumanMessage, AIMessage

# --- Configuration ---
# MODEL_NAME = "gemini-3-flash-preview"
# MODEL_PROVIDER = "google_genai"
MODEL_NAME = "my-model"
MODEL_PROVIDER = "openai"
BASE_URL = "https://myung-courageous-gino.ngrok-free.dev/v1"
API_KEY = "not-needed"

llm = init_chat_model(
    model=MODEL_NAME,
    model_provider=MODEL_PROVIDER,
    base_url=BASE_URL,
    api_key=API_KEY,
    temperature=0.2,
)

# Step 1: Create the base prompt template
base_prompt = ChatPromptTemplate.from_messages([
    ("system",
     "You are a medical specialist in {specialty}. "
     "Answer patient questions from the perspective of your specialty. "
     "Keep answers concise (under 150 words)."),
    MessagesPlaceholder("history"),
    ("human", "{question}"),
])

# Step 2: Create specialist prompts using partial()
cardiology_prompt = base_prompt.partial(specialty="cardiology")
neurology_prompt = base_prompt.partial(specialty="neurology")


# Step 3: make_specialist_chain() function
def make_specialist_chain(specialty: str):
    """Create a chain for a given medical specialty."""
    return base_prompt.partial(specialty=specialty) | llm


# Step 4: Ask the same question to two specialists (2 API calls)
patient_question = (
    "I have been experiencing chest pain and shortness of breath. "
    "What could be causing this?"
)

specialties = ["cardiology", "neurology"]
responses = {}  # save responses so we can reuse cardiology's for history

for spec in specialties:
    chain = make_specialist_chain(spec)
    response = chain.invoke({"history": [], "question": patient_question})
    responses[spec] = response

    print(f"--- {spec.upper()} specialist ---")
    print(response.text)
    print(f"Tokens used: {response.usage_metadata}")
    print()

# Step 5: Follow-up with cardiologist using conversation history (1 API call)
print("=" * 60)
print("FOLLOW-UP WITH CARDIOLOGIST")
print("=" * 60)

cardiology_chain = make_specialist_chain("cardiology")

# Build history from the first cardiology exchange
history = [
    HumanMessage(content=patient_question),
    responses["cardiology"],
]

question_followup = "What tests would you recommend?"
response_followup = cardiology_chain.invoke(
    {"history": history, "question": question_followup}
)
print(f"\nPatient: {question_followup}")
print(f"Cardiologist: {response_followup.text}")
print(f"Tokens: {response_followup.usage_metadata}")
