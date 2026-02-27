"""
Solution — MCP Math Server (extended)
=======================================
MCP server with add, multiply, subtract, and power tools.

This file is used by the MCP client. Do NOT run it directly.
"""

from mcp.server.fastmcp import FastMCP

mcp = FastMCP("Math")

@mcp.tool()
def add(a: int, b: int) -> int:
    """Add two numbers"""
    return a + b

@mcp.tool()
def multiply(a: int, b: int) -> int:
    """Multiply two numbers"""
    return a * b

@mcp.tool()
def subtract(a: int, b: int) -> int:
    """Subtract the second number from the first"""
    return a - b

@mcp.tool()
def power(a: int, b: int) -> int:
    """Raise the first number to the power of the second"""
    return a ** b

if __name__ == "__main__":
    mcp.run(transport="stdio")
