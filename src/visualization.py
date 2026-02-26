import os
import sys
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import config

def plot_alignment_scores(score_dict, filename="alignment_scores.png"):
    """
    Creates a bar chart comparing different alignment scores.
    """
    labels = list(score_dict.keys())
    scores = list(score_dict.values())
    
    plt.figure(figsize=(8, 5))
    plt.bar(labels, scores, color=['#4CAF50', '#2196F3', '#FFC107'])
    plt.title('Comparison of Sequence Alignment Scores')
    plt.ylabel('Score')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    output_path = os.path.join(config.FIGURES_DIR, filename)
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved successfully to {output_path}")
    plt.close()

if __name__ == "__main__":
    # Dummy data based on typical NW vs SW score variations
    sample_scores = {
        "Global (Needleman-Wunsch)": 150.5,
        "Local (Smith-Waterman)": 165.0,
        "ClustalW SP Score": 142.0
    }
    plot_alignment_scores(sample_scores)