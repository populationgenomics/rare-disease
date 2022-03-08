# Design Doc

## Aim

Take in variant data, and run a variant prioritisation algorithm. Highlight a minimal set of variants per sample, based on strict criteria

In general this is focused on True Positives, and is happy to allow False Negatives, making a highly specific actionable dataset, at the potential expense of sensitivity

This is currently in an MVP phase, where the core logic should be accurate but the input and output are only rough drafts

## Concept

Take in variant data