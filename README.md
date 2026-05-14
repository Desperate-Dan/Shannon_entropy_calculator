# Shannon_entropy_calculator
Python script to calculate and plot the Shannon entropy for a DNA alignment.

## Installation
```bash
pip install plotly
```

## Usage
```
python entropy_calculator.py <alignment.fasta> [title] [smooth_factor] [--lower-threshold VALUE] [--upper-threshold VALUE]
```

## Parameters
- `<alignment.fasta>` (required): Input FASTA alignment file (supports .gz)
- `[title]` (optional): Output filename for the graph. Defaults to alignment filename. Output will be in HTML format
- `[smooth_factor]` (optional): Smoothing window size. Defaults to 0 (no smoothing). A value of 10 is a good starting place
- `[--lower-threshold VALUE]` (optional): Highlight regions where entropy is **below or equal to** this value (red shading)
- `[--upper-threshold VALUE]` (optional): Highlight regions where entropy is **above or equal to** this value (red shading)

## Examples
```bash
# Basic analysis (outputs alignment.html)
python entropy_calculator.py alignment.fasta

# With smoothing
python entropy_calculator.py alignment.fasta output.html 10

# Highlight low entropy regions (< 0.5)
python entropy_calculator.py alignment.fasta output.html 0 --lower-threshold 0.5

# Highlight high entropy regions (> 2.0)
python entropy_calculator.py alignment.fasta output.html 0 --upper-threshold 2.0

# Combined: smooth + highlight
python entropy_calculator.py alignment.fasta output.html 10 --lower-threshold 1.0
```

## Features
- Calculates Shannon entropy for each position in a DNA alignment
- Generates **interactive HTML plots** with full zoom and pan capabilities
- Range slider at the bottom for quick navigation across the alignment
- Hover over points to see exact position and entropy values
- Displays max entropy reference line (2.321 for 5 bases)
- Highlights threshold regions with red shading for easy identification of conservation/variation
- Responsive design with high-quality output

