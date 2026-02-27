"""
Solution — Exercise 1: Ollama with LangChain
"""

from langchain.chat_models import init_chat_model
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser

# Step 1: Create the model
llm = init_chat_model(
    model="qwen3:0.6b",
    model_provider="ollama",
    temperature=0.2,
)

# Step 2: Create the prompt template
prompt = ChatPromptTemplate.from_messages([
    ("system", "You are a helpful biology tutor. Explain topics clearly and "
               "in simple terms that any student can understand."),
    ("human", "{question}"),
])

# Step 3: Build the chain
chain = prompt | llm | StrOutputParser()

# Step 4: Invoke and print
result = chain.invoke({"question": "What is photosynthesis?"})
print(result)
