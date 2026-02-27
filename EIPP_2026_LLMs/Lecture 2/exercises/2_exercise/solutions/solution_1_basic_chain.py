"""
Solution — Exercise 1: Basic Chain Construction
"""

from dotenv import load_dotenv
load_dotenv()

from langchain.chat_models import init_chat_model
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser

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

# Step 1: Create the prompt template
prompt = ChatPromptTemplate.from_messages([
    ("system", "You are a friendly science teacher. Explain complex scientific "
               "concepts in simple terms that a high school student can understand."),
    ("human", "{concept}"),
])

# Step 2: Build the chain
chain = prompt | llm | StrOutputParser()

# Step 3: Invoke and print
result = chain.invoke({"concept": "photosynthesis"})
print(result)
