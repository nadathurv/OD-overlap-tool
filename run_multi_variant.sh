#!/usr/bin/env bash
set -e

echo "🔄 Running multi-variant CDSCO vs FDA comparison..."
echo "This will create separate overlap files for each CDSCO variant:"
echo "  - overlap_cdsco_clean_names_only.csv"
echo "  - overlap_unclean_cdsco.csv"
echo "  - overlap_vikram_cdsco_clean.csv"
echo ""

# Run unified script in multi-variant mode
python compare.py --multi-variant "$@"

echo ""
echo "✅ Multi-variant comparison complete!"
echo "📁 Results saved in data/processed/"
echo ""
echo "Summary of matches found:"
if [ -f "data/processed/overlap_cdsco_clean_names_only.csv" ]; then
    names_only_count=$(tail -n +2 "data/processed/overlap_cdsco_clean_names_only.csv" | wc -l)
    echo "  - Names only variant: $names_only_count matches"
fi

if [ -f "data/processed/overlap_unclean_cdsco.csv" ]; then
    unclean_count=$(tail -n +2 "data/processed/overlap_unclean_cdsco.csv" | wc -l)
    echo "  - Unclean variant: $unclean_count matches"
fi

if [ -f "data/processed/overlap_vikram_cdsco_clean.csv" ]; then
    vikram_count=$(tail -n +2 "data/processed/overlap_vikram_cdsco_clean.csv" | wc -l)
    echo "  - Vikram cleaned variant: $vikram_count matches"
fi 