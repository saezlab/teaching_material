"""
MCP Math Server
================
A minimal MCP server exposing math tools over stdio.
Students extend this with subtract and power tools.

This file is used by the MCP client exercise.
Do NOT run this file directly — the client starts it automatically.
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

# TODO: Add a subtract(a: int, b: int) tool that returns a - b
# Don't forget the @mcp.tool() decorator and a docstring!


# TODO: Add a power(a: int, b: int) tool that returns a ** b
# Don't forget the @mcp.tool() decorator and a docstring!


if __name__ == "__main__":
    mcp.run(transport="stdio")
