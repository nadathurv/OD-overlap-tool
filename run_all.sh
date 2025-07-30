#!/usr/bin/env bash
set -e

echo "🧬 Orphan Drug Overlap Pipeline"
echo "================================"

# Check for explode flag
EXPLODE_FLAG=""
if [[ "$1" == "--explode-combinations" ]]; then
    EXPLODE_FLAG="--use-exploded"
    echo "✓ Using exploded combination drugs"
fi

# Step 1: Clean raw data
echo ""
echo "📊 Step 1: Cleaning raw data..."
python -m src.data.clean $EXPLODE_FLAG

# Step 2: Run comparison
echo ""
echo "🔍 Step 2: Running drug comparison..."
python compare.py $EXPLODE_FLAG

# Summary
echo ""
echo "✅ Pipeline complete!"
echo "📁 Results saved to: data/processed/overlap.csv"

# Show match count
if [ -f "data/processed/overlap.csv" ]; then
    MATCHES=$(tail -n +2 "data/processed/overlap.csv" | wc -l)
    echo "🎯 Total matches found: $MATCHES"
fi
