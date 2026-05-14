import math
import sys
import gzip
import re
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def shannon_entropy(column):
    """Calculate Shannon entropy for a column of nucleotides"""
    counts = {}
    total = 0
    for base in column:
        """I don't think N should be in this list as it already means any base so is uninformative"""
        if base not in 'ACGT-':
            continue  # skip invalid characters
        counts[base] = counts.get(base, 0) + 1
        total += 1
    if total == 0:
        return 0.0
    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)

    return entropy

def read_alignment(filename):
    """Check if alignment file is gzipped or not"""
    if re.search(".gz", filename):
        with gzip.open(filename,"rt") as f:
            sequences = sequence_reader(f)
            return(sequences)
    else:
        with open(filename,"r") as f:
            sequences = sequence_reader(f)
            return(sequences)
    

def sequence_reader(f,sequences=[]):
    """Read alignment file in FASTA format and return list of sequences"""   
    seq = ''
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if seq:
                sequences.append(seq)
                seq = ''
        else:
            seq += line.upper()
    if seq:
        sequences.append(seq)
    return sequences

def calculate_entropy_matrix(sequences):
    """Calculate entropy for each column in the alignment"""
    if not sequences:
        return []
    length = len(sequences[0])
    entropies = []
    for i in range(length):
        column = [seq[i] for seq in sequences if len(seq) > i]
        entropies.append(shannon_entropy(column))
    
    return entropies

def plotter(entropies, title, smooth, lower_threshold=None, upper_threshold=None):
    """Plotting the entropies in an interactive line graph with optional threshold highlighting"""
    entropies_dict = {}
    for n, i in enumerate(entropies):
        if not smooth:
            entropies_dict[n] = i
        else:
            """Smoothing works by averageing the entropy across a window of size smooth"""
            if n < smooth:
                entropies_dict[n] = (sum(entropies[0:(n+int(smooth/2))]))/(n+int(smooth/2))
            elif n > (len(entropies) - smooth):
                entropies_dict[n] = (sum(entropies[(n-int(smooth/2)):len(entropies)]))/((len(entropies)-n)+(int(smooth/2)))
            else:
                if smooth % 2 == 0:
                    entropies_dict[n] = (sum(entropies[(n-int(smooth/2)):(n+int(smooth/2))]))/(smooth)
                else:
                    entropies_dict[n] = (sum(entropies[(n-int(smooth/2)):(1+n+int(smooth/2))]))/(smooth)
    
    fig = go.Figure()
    
    # Add entropy line
    positions = sorted(entropies_dict.keys())
    values = [entropies_dict[p] for p in positions]
    
    fig.add_trace(go.Scatter(
        x=positions,
        y=values,
        mode='lines',
        name='Shannon Entropy',
        line=dict(color='blue', width=2),
        hovertemplate='Position: %{x}<br>Entropy: %{y:.3f}<extra></extra>'
    ))
    
    # Add max entropy reference line
    fig.add_hline(
        y=2.321,
        line_dash="dash",
        line_color="gray",
        annotation_text="Max entropy (5 bases)",
        annotation_position="right"
    )
    
    # Highlight regions based on thresholds
    if lower_threshold is not None or upper_threshold is not None:
        add_threshold_regions(fig, entropies_dict, lower_threshold, upper_threshold)
    
    # Update layout for interactivity
    fig.update_layout(
        title=title,
        xaxis_title="Amplified Region (Position)",
        yaxis_title="Shannon Entropy",
        hovermode='x unified',
        width=1200,
        height=600,
        template='plotly_white',
        yaxis=dict(range=[0, 2.5]),
        xaxis=dict(
            rangeslider=dict(visible=True, thickness=0.05),
            type="linear"
        )
    )
    
    # Save as interactive HTML
    if not title.endswith('.html'):
        title = title.replace('.png', '') + '.html'
    
    fig.write_html(title)
    print(f"Interactive plot saved to: {title}")

def add_threshold_regions(fig, entropies_dict, lower_threshold=None, upper_threshold=None):
    """Add highlighted regions where entropy meets threshold criteria"""
    positions = sorted(entropies_dict.keys())
    if not positions:
        return
    
    # Find continuous regions that meet threshold criteria
    regions = []
    in_region = False
    region_start = None
    
    for pos in positions:
        entropy_val = entropies_dict[pos]
        meets_criteria = False
        
        if lower_threshold is not None and upper_threshold is not None:
            meets_criteria = lower_threshold <= entropy_val <= upper_threshold
        elif lower_threshold is not None:
            meets_criteria = entropy_val <= lower_threshold
        elif upper_threshold is not None:
            meets_criteria = entropy_val >= upper_threshold
        
        if meets_criteria:
            if not in_region:
                region_start = pos
                in_region = True
        else:
            if in_region:
                regions.append((region_start, pos - 1))
                in_region = False
    
    # Don't forget the last region if it extends to the end
    if in_region:
        regions.append((region_start, positions[-1]))
    
    # Add shaded regions
    for i, (start, end) in enumerate(regions):
        fig.add_vrect(
            x0=start, x1=end,
            fillcolor="red", opacity=0.15,
            layer="below", line_width=0
        )

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python entropy_calculator.py <alignment.fasta> [title] [smooth_factor] [--lower-threshold VALUE] [--upper-threshold VALUE]")
        sys.exit(1)
    sequences = read_alignment(sys.argv[1])
    try: 
        title = sys.argv[2]
    except:
        title = sys.argv[1].split(".")[0]
    try: 
        smooth = int(sys.argv[3])
    except:
        smooth = 0
    
    # Parse threshold arguments
    lower_threshold = None
    upper_threshold = None
    try:
        if '--lower-threshold' in sys.argv:
            lower_idx = sys.argv.index('--lower-threshold')
            lower_threshold = float(sys.argv[lower_idx + 1])
        if '--upper-threshold' in sys.argv:
            upper_idx = sys.argv.index('--upper-threshold')
            upper_threshold = float(sys.argv[upper_idx + 1])
    except (ValueError, IndexError):
        print("Error: threshold values must be valid numbers")
        sys.exit(1)
    
    entropies = calculate_entropy_matrix(sequences)
    plotter(entropies, title, smooth, lower_threshold, upper_threshold)
