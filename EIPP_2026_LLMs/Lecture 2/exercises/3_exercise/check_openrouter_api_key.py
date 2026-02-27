"""Verify OpenRouter API key is set correctly."""
import os
from dotenv import load_dotenv

load_dotenv()

if "OPENROUTER_API_KEY" not in os.environ:
    print("❌ OPENROUTER_API_KEY not found in environment")
    print("   Make sure you added it to playground/.env")
    exit(1)

key = os.environ["OPENROUTER_API_KEY"]
if not key.startswith("sk-or-"):
    print("⚠️  Warning: OpenRouter keys typically start with 'sk-or-'")
    print(f"   Your key starts with: {key[:10]}")

print("✅ OPENROUTER_API_KEY is set")
print(f"   Key prefix: {key[:15]}...")
