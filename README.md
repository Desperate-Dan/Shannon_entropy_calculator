# Shannon_entropy_calculator
Python script to calculate and plot the Shannon entropy for a DNA alignment.

## Usage
```python entropy_calculator.py <alignment.fasta> <title> <smooth_factor>```

These are positional arguments so you need a title if you're using the smoothing factor. Smoothing averages the entropy values across a window of size 'smooth'. Off by default, but a value of 10 is a good starting place.
