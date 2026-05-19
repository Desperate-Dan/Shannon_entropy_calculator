import streamlit as st
import plotly.graph_objects as go
from Shannon_entropy_calculator import read_alignment, calculate_entropy_matrix
import os

# Page configuration
st.set_page_config(page_title="Shannon Entropy Calculator", layout="wide", initial_sidebar_state="expanded")

st.title("📊 Shannon Entropy Calculator - Interactive Analysis")

# Sidebar for file selection and controls
with st.sidebar:
    st.header("Controls")
    
    # File selection
    st.subheader("1. Select Alignment File")
    
    # List available FASTA files
    fasta_files = [f for f in os.listdir('.') if f.endswith(('.fasta', '.fa', '.aln', '.fasta.gz', '.fa.gz'))]
    
    if fasta_files:
        selected_file = st.selectbox("Choose from existing files:", fasta_files)
    else:
        selected_file = None
        st.warning("No FASTA files found in current directory")
    
    # File upload option
    uploaded_file = st.file_uploader("Or upload a new alignment file:", type=['fasta', 'fa', 'aln', 'fasta.gz', 'fa.gz'])
    
    if uploaded_file is not None:
        # Save uploaded file temporarily
        with open(uploaded_file.name, 'wb') as f:
            f.write(uploaded_file.getbuffer())
        alignment_file = uploaded_file.name
        st.success(f"Uploaded: {uploaded_file.name}")
    elif selected_file:
        alignment_file = selected_file
    else:
        alignment_file = None
        st.error("Please upload or select an alignment file")

# Main content
if alignment_file:
    try:
        # Load the alignment
        sequences = read_alignment(alignment_file)
        entropies = calculate_entropy_matrix(sequences)
        
        # Display file info
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Alignment Length", len(entropies))
        with col2:
            st.metric("Number of Sequences", len(sequences))
        with col3:
            st.metric("Max Entropy", f"{max(entropies):.3f}")
        
        st.divider()
        
        # Parameter controls in sidebar
        with st.sidebar:
            st.subheader("2. Adjust Parameters")
            
            # Smoothing factor
            smooth = st.slider(
                "Smoothing Factor",
                min_value=0,
                max_value=min(100, len(entropies) // 10),
                value=5,
                step=1,
                help="Window size for smoothing. 0 = no smoothing"
            )
            
            st.subheader("3. Threshold Highlighting")
            
            # Threshold options
            threshold_type = st.radio(
                "Highlight regions:",
                ["None", "Below threshold", "Above threshold", "Between thresholds"]
            )
            
            lower_threshold = None
            upper_threshold = None
            
            if threshold_type == "Below threshold":
                lower_threshold = st.slider(
                    "Lower threshold",
                    min_value=0.0,
                    max_value=2.5,
                    value=1.0,
                    step=0.1,
                    help="Highlight regions with entropy ≤ this value"
                )
            elif threshold_type == "Above threshold":
                upper_threshold = st.slider(
                    "Upper threshold",
                    min_value=0.0,
                    max_value=2.5,
                    value=2.0,
                    step=0.1,
                    help="Highlight regions with entropy ≥ this value"
                )
            elif threshold_type == "Between thresholds":
                col1, col2 = st.columns(2)
                with col1:
                    lower_threshold = st.slider(
                        "Lower threshold",
                        min_value=0.0,
                        max_value=2.5,
                        value=1.0,
                        step=0.1,
                        key="lower"
                    )
                with col2:
                    upper_threshold = st.slider(
                        "Upper threshold",
                        min_value=0.0,
                        max_value=2.5,
                        value=2.0,
                        step=0.1,
                        key="upper"
                    )
        
        # Calculate smoothed entropies
        entropies_dict = {}
        for n, i in enumerate(entropies):
            if not smooth:
                entropies_dict[n] = i
            else:
                if n < smooth:
                    entropies_dict[n] = (sum(entropies[0:(n+int(smooth/2))]))/(1+n+int(smooth/2))
                elif n > (len(entropies) - smooth):
                    entropies_dict[n] = (sum(entropies[(n-int(smooth/2)):len(entropies)]))/((len(entropies)-n)+(int(smooth/2)))
                else:
                    entropies_dict[n] = (sum(entropies[(n-int(smooth/2)):(n+int(smooth/2))]))/(smooth)
        
        # Create interactive plot
        fig = go.Figure()
        
        # Add entropy line
        positions = sorted(entropies_dict.keys())
        values = [entropies_dict[p] for p in positions]
        
        fig.add_trace(go.Scatter(
            x=positions,
            y=values,
            mode='lines',
            name='Shannon Entropy',
            line=dict(color='#1f77b4', width=2),
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
        
        # Highlight threshold regions if selected
        if threshold_type != "None":
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
            
            if in_region:
                regions.append((region_start, positions[-1]))
            
            # Add shaded regions
            for start, end in regions:
                fig.add_vrect(
                    x0=start, x1=end,
                    fillcolor="red", opacity=0.15,
                    layer="below", line_width=0
                )
            
            # Show statistics about highlighted regions
            total_highlighted = sum(end - start + 1 for start, end in regions)
            st.sidebar.info(f"📍 {len(regions)} region(s) highlighted\n{total_highlighted} positions meet criteria")
        
        # Update layout
        fig.update_layout(
            title=f"Shannon Entropy Analysis - {os.path.basename(alignment_file)}",
            xaxis_title="Position in Alignment",
            yaxis_title="Shannon Entropy",
            hovermode='x unified',
            height=600,
            template='plotly_white',
            yaxis=dict(range=[0, 2.5]),
            xaxis=dict(rangeslider=dict(visible=True, thickness=0.05))
        )
        
        # Display plot
        st.plotly_chart(fig, use_container_width=True)
        
        # Export options
        st.divider()
        with st.expander("📥 Export Options"):
            col1, col2 = st.columns(2)
            
            with col1:
                if st.button("Save as HTML"):
                    output_name = f"{os.path.splitext(alignment_file)[0]}_analysis.html"
                    fig.write_html(output_name)
                    st.success(f"Saved to {output_name}")
            
            with col2:
                if st.button("Save as PNG"):
                    try:
                        output_name = f"{os.path.splitext(alignment_file)[0]}_analysis.png"
                        fig.write_image(output_name)
                        st.success(f"Saved to {output_name}")
                    except:
                        st.warning("PNG export requires kaleido. Install with: pip install kaleido")
            
            # Download entropy values as CSV
            st.subheader("Download Data")
            csv_data = "Position,Entropy\n"
            for pos in positions:
                csv_data += f"{pos},{entropies_dict[pos]:.6f}\n"
            
            st.download_button(
                label="Download entropy values (CSV)",
                data=csv_data,
                file_name=f"{os.path.splitext(alignment_file)[0]}_entropy.csv",
                mime="text/csv"
            )
            
            # Export highlighted regions if thresholds are applied
            if threshold_type != "None":
                st.subheader("Export Highlighted Regions")
                
                # Identify highlighted regions
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
                
                if in_region:
                    regions.append((region_start, positions[-1]))
                
                # Generate CSV with expanded regions
                region_csv = "Region,Position,Entropy\n"
                for region_num, (start, end) in enumerate(regions, 1):
                    # Expand region by smooth factor
                    expanded_start = max(0, start - int(smooth/2))
                    expanded_end = min(positions[-1], end + int(smooth/2))
                    
                    for pos in range(expanded_start, expanded_end + 1):
                        if pos in entropies_dict:
                            region_csv += f"region_{region_num},{pos},{entropies_dict[pos]:.6f}\n"
                
                st.download_button(
                    label="Download highlighted regions (CSV)",
                    data=region_csv,
                    file_name=f"{os.path.splitext(alignment_file)[0]}_highlighted_regions.csv",
                    mime="text/csv"
                )
        
        # Summary statistics
        with st.expander("📈 Summary Statistics"):
            col1, col2, col3, col4, col5 = st.columns(5)
            
            import statistics
            values_list = list(values)
            
            with col1:
                st.metric("Mean", f"{statistics.mean(values_list):.3f}")
            with col2:
                st.metric("Median", f"{statistics.median(values_list):.3f}")
            with col3:
                st.metric("Std Dev", f"{statistics.stdev(values_list):.3f}")
            with col4:
                st.metric("Min", f"{min(values_list):.3f}")
            with col5:
                st.metric("Max", f"{max(values_list):.3f}")
    
    except Exception as e:
        st.error(f"Error processing file: {str(e)}")
        st.info("Make sure the file is a valid FASTA alignment")

else:
    st.info("👈 Please select or upload an alignment file to begin analysis")
