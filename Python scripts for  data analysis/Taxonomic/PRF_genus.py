import pandas as pd
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate precision, recall, and F-score for various genera from a given dataset.')
    parser.add_argument('file_path', type=str, help='The path to the CSV file to be processed.')
    return parser.parse_args()

def calculate_metrics(df):
    results = {}
    genera = ["LF", "BS", "SA", "PA", "SE", "EC", "EF", "LM", "CN", "SC"]
    for genus in genera:
        tp = df.loc[df["Category"] == f"D_{genus}", "Count"].values[0]
        fp = df.loc[df["Category"] == f"B_{genus}", "Count"].values[0] - df.loc[df["Category"] == f"D_{genus}", "Count"].values[0]
        fn = df.loc[df["Category"] == f"A_{genus}", "Count"].values[0] - df.loc[df["Category"] == f"D_{genus}", "Count"].values[0]
        
        precision = tp / (tp + fp) if tp + fp > 0 else 0
        recall = tp / (tp + fn) if tp + fn > 0 else 0
        f_score = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0
        
        results[genus] = {
            "Precision": precision,
            "Recall": recall,
            "F-Score": f_score
        }
    return results

def main():
    args = parse_arguments()
    df = pd.read_csv(args.file_path)  # Assumes the file is in CSV format

    metrics = calculate_metrics(df)
    
    # Print results
    for genus, metrics in metrics.items():
        print(f"{genus}: Precision = {metrics['Precision']:.4f}, Recall = {metrics['Recall']:.4f}, F-Score = {metrics['F-Score']:.4f}")

if __name__ == "__main__":
    main()
