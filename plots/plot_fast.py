#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def main():
    # Load the CSV file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    csv_path = os.path.join(project_dir, 'feeddown_dispose/Proton/finalise_feeddown_dispose_default.csv')
    
    # Read the data
    df = pd.read_csv(csv_path)
    
    # Filter data for Intg differential type and SS/OS pair types
    intg_data = df[(df['diff_type'] == 'Intg') & (df['diff_bin'] == 0.5) & 
                   ((df['pair_type'] == 'SS') | (df['pair_type'] == 'OS'))]
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Define colors and markers
    colors = {'SS': 'blue', 'OS': 'red'}
    markers = {'SS': 'o', 'OS': 's'}
    
    # Group by pair_type
    for pair_type, group in intg_data.groupby('pair_type'):
        # Sort by centrality
        group = group.sort_values('centrality')
        
        # Plot delta vs centrality
        ax1.errorbar(
            group['centrality'], 
            group['delta'], 
            yerr=group['delta_err'], 
            fmt=markers[pair_type],
            color=colors[pair_type],
            label=f'{pair_type}',
            capsize=3,
            markersize=6,
            elinewidth=1,
            linewidth=1.5
        )
        
        # Plot gamma vs centrality
        ax2.errorbar(
            group['centrality'], 
            group['gamma'], 
            yerr=group['gamma_err'], 
            fmt=markers[pair_type],
            color=colors[pair_type],
            label=f'{pair_type}',
            capsize=3,
            markersize=6,
            elinewidth=1,
            linewidth=1.5
        )
    
    # Set titles and labels
    ax1.set_title(r'$\delta$ vs Centrality for Intg')
    ax1.set_xlabel('Centrality (%)')
    ax1.set_ylabel(r'$\delta$')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    ax2.set_title(r'$\gamma$ vs Centrality for Intg')
    ax2.set_xlabel('Centrality (%)')
    ax2.set_ylabel(r'$\gamma$')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Add a main title
    plt.suptitle('Delta and Gamma vs Centrality for Integrated Differential Type', fontsize=14)
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    
    # Save the figure
    output_dir = os.path.join(project_dir, 'plots/output')
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, 'intg_delta_gamma_vs_centrality.png'), dpi=300)
    
    # Also save as PDF for publications
    plt.savefig(os.path.join(output_dir, 'intg_delta_gamma_vs_centrality.pdf'))
    
    print(f"Plots saved to {output_dir}")
    
    # Display the plot
    plt.show()

if __name__ == "__main__":
    main()