"""
Solution — Exercise 2: Message Types and Conversation
"""

from dotenv import load_dotenv
load_dotenv()

from langchain.chat_models import init_chat_model
from langchain_core.messages import SystemMessage, HumanMessage

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

# Step 1: Create the system message
system_msg = SystemMessage(
    content="You are a knowledgeable and empathetic doctor. Answer patient "
            "questions about medications clearly and reassuringly. Keep your "
            "answers concise (under 150 words)."
)

# Step 2: First exchange
messages = [system_msg]

patient_q1 = HumanMessage(content="What are the common side effects of ibuprofen?")
messages.append(patient_q1)

response_1 = llm.invoke(messages)
messages.append(response_1)

print("--- Doctor's first response ---")
print(response_1.text)

# Step 3: Second exchange (follow-up)
patient_q2 = HumanMessage(
    content="Should I be worried if I take it daily for a week?"
)
messages.append(patient_q2)

response_2 = llm.invoke(messages)
messages.append(response_2)

print("\n--- Doctor's second response ---")
print(response_2.text)

# Step 4: Print message types
print("\n--- Message types in conversation ---")
for msg in messages:
    print(type(msg))
