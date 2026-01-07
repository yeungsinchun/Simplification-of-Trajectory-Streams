import pandas as pd
import sys

def analyze_file(file_path, metric_col, baseline_col, reduction_col_name, better_condition='smaller'):
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return None

    # Group by ID and pick the best row based on the metric
    # For points, we want min points. For dist, we want min dist.
    if better_condition == 'smaller':
        best_rows = df.loc[df.groupby('id')[metric_col].idxmin()]
    else:
        # Just in case we ever want larger, but for now we want smaller for both
        best_rows = df.loc[df.groupby('id')[metric_col].idxmax()]

    # Calculate reduction percentage for the selected best rows
    # (Baseline - New) / Baseline * 100
    best_rows[reduction_col_name] = (best_rows[baseline_col] - best_rows[metric_col]) / best_rows[baseline_col] * 100
    
    return best_rows

def main():
    # 1. Analyze Size Reduction (minimize_points.csv)
    print("--- Analysis of Size Reduction (minimize_points.csv) ---")
    points_df = analyze_file("minimize_points.csv", "best_simp_points", "dp_points", "points_reduction_pct")
    
    if points_df is not None:
        avg_dp_points = points_df['dp_points'].mean()
        avg_simp_points = points_df['best_simp_points'].mean()
        
        better_cases = points_df[points_df['best_simp_points'] <= points_df['dp_points']]
        num_better = len(better_cases)
        total = len(points_df)
        pct_better = (num_better / total) * 100
        
        avg_reduction = better_cases['points_reduction_pct'].mean() if num_better > 0 else 0

        print(f"Total Trajectories: {total}")
        print(f"Average DP Points: {avg_dp_points:.2f}")
        print(f"Average Simp Points (Best): {avg_simp_points:.2f}")
        print(f"Cases with Reduced Size: {num_better}/{total} ({pct_better:.1f}%)")
        print(f"Average Size Reduction (where reduced): {avg_reduction:.2f}%")
        print("\n")

    # 2. Analyze Error Reduction (minimize_frechet.csv)
    print("--- Analysis of Error Reduction (minimize_frechet.csv) ---")
    frechet_df = None
    valid_id_pct = 0
    valid_ids_count = 0
    total_ids = 0
    
    try:
        raw_df = pd.read_csv("minimize_frechet.csv")
        total_ids = raw_df['id'].nunique()
        
        # Filter: only consider rows where simp_points <= dp_points
        valid_rows = raw_df[raw_df['simp_points'] <= raw_df['dp_points']]
        
        valid_ids_count = valid_rows['id'].nunique()
        valid_id_pct = (valid_ids_count / total_ids) * 100 if total_ids > 0 else 0
        
        print(f"Total Trajectories: {total_ids}")
        print(f"Trajectories with Simp Points <= DP Points: {valid_ids_count}/{total_ids} ({valid_id_pct:.1f}%)")
        
        if valid_ids_count > 0:
            # For each ID, pick the row with the minimum best_simp_dist among valid rows
            frechet_df = valid_rows.loc[valid_rows.groupby('id')['best_simp_dist'].idxmin()].copy()
            
            # Calculate reductions
            frechet_df['dist_reduction_pct'] = (frechet_df['dp_dist'] - frechet_df['best_simp_dist']) / frechet_df['dp_dist'] * 100
            frechet_df['points_reduction_pct'] = (frechet_df['dp_points'] - frechet_df['simp_points']) / frechet_df['dp_points'] * 100
            
            avg_dist_reduction = frechet_df['dist_reduction_pct'].mean()
            avg_points_reduction = frechet_df['points_reduction_pct'].mean()
            
            avg_dp_dist = frechet_df['dp_dist'].mean()
            avg_simp_dist = frechet_df['best_simp_dist'].mean()
            avg_dp_points = frechet_df['dp_points'].mean()
            avg_simp_points = frechet_df['simp_points'].mean()
            
            print(f"Average DP Distance: {avg_dp_dist:.2f}")
            print(f"Average Simp Distance (Best Valid): {avg_simp_dist:.2f}")
            print(f"Average Error Reduction (for valid cases): {avg_dist_reduction:.2f}%")
            
            print(f"Average DP Points: {avg_dp_points:.2f}")
            print(f"Average Simp Points (Best Valid): {avg_simp_points:.2f}")
            print(f"Average Points Reduction (for valid cases): {avg_points_reduction:.2f}%")
            
            # Check for anomalies
            worse_cases = frechet_df[frechet_df['best_simp_dist'] > frechet_df['dp_dist']]
            num_worse = len(worse_cases)
            if num_worse > 0:
                print(f"\nWARNING: Found {num_worse} valid cases where Simp Dist > DP Dist:")
                print(worse_cases[['id', 'dp_dist', 'best_simp_dist', 'dp_points', 'simp_points']].to_string(index=False))
                
    except FileNotFoundError:
        print("Error: File 'minimize_frechet.csv' not found.")

    with open("analysis_summary.txt", "w") as f:
        if points_df is not None:
            f.write("--- Size Reduction ---\n")
            f.write(f"Cases with Reduced Size: {len(points_df[points_df['best_simp_points'] < points_df['dp_points']])}/{len(points_df)} ({len(points_df[points_df['best_simp_points'] < points_df['dp_points']])/len(points_df)*100:.1f}%)\n")
            f.write(f"Average Size Reduction: {points_df[points_df['best_simp_points'] < points_df['dp_points']]['points_reduction_pct'].mean():.2f}%\n\n")
        
        if frechet_df is not None:
            f.write("--- Error Reduction (Valid Cases: Simp Points <= DP Points) ---\n")
            f.write(f"Trajectories with Simp Points <= DP Points: {valid_ids_count}/{total_ids} ({valid_id_pct:.1f}%)\n")
            f.write(f"Average Error Reduction: {frechet_df['dist_reduction_pct'].mean():.2f}%\n")
            f.write(f"Average Points Reduction: {frechet_df['points_reduction_pct'].mean():.2f}%\n")
            f.write(f"Average DP Points: {frechet_df['dp_points'].mean():.2f}\n")
            f.write(f"Average Simp Points: {frechet_df['simp_points'].mean():.2f}\n")

if __name__ == "__main__":
    main()
