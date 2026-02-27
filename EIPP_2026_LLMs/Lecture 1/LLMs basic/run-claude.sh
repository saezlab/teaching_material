#!/bin/bash -l
set -e

PROMPT_FILE="$1"
ITERATIONS="${2:-0}"  # 0 means infinite

if [ ! -f "$PROMPT_FILE" ]; then
    echo "Error: Prompt file '$PROMPT_FILE' not found"
    exit 1
fi

run_iteration() {
    local i=$1
    echo ""
    echo "========================================"
    echo "Iteration $i - $(date '+%Y-%m-%d %H:%M:%S')"
    echo "========================================"

    cat "$PROMPT_FILE" | claude -p \
        --dangerously-skip-permissions \
        --output-format=stream-json \
        --model opus \
        --verbose \
    | jq -r '
        if .type == "assistant" then
            .message.content[] |
            if .type == "text" then
                "💬 \(.text)"
            elif .type == "tool_use" then
                "🔧 \(.name): \(.input | tostring | .[0:200])"
            else empty end
        elif .type == "result" then
            "\n📊 Done | Cost: $\(.total_cost_usd | tostring | .[0:8]) | Duration: \(.duration_ms)ms | Turns: \(.num_turns)"
        else empty end
    '
}

i=1
while true; do
    run_iteration $i

    if [ "$ITERATIONS" -gt 0 ] && [ "$i" -ge "$ITERATIONS" ]; then
        echo ""
        echo "Completed $ITERATIONS iterations"
        break
    fi

    ((i++))
done
