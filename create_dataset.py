import os
import re
import pandas as pd

def analyze_files(base_folder, islets=None, glucose_levels=None, capillaries=None):
    """
    Analyze files in folders based on user-defined criteria.

    Parameters:
        base_folder (str): The base folder containing subfolders with data files.
        islets (list of str): List of islet identifiers to filter (e.g., ['H51', 'H52']).
        glucose_levels (list of str): List of glucose levels to filter (e.g., ['G6', 'G7']).
        capillaries (list of int): List of capillary counts to filter (e.g., [5, 10]).

    Returns:
        pd.DataFrame: A DataFrame containing the analysis results.
    """
    results = []

    # Walk through the base folder and its subfolders
    for root, dirs, files in os.walk(base_folder):
        for file in files:
            if file.endswith('_data_test.dat') and file.startswith('H5'):
                file_path = os.path.join(root, file)

                # Extract metadata from the folder structure or filename
                folder_name = os.path.basename(root)
                file_name = os.path.splitext(file)[0]

                # Use regex to extract islet, glucose level, and capillary count
                islet_match = re.search(r'(H\d+)', folder_name)
                glucose_match = re.search(r'(G\d+)', folder_name)
                capillary_match = re.search(r'(\d+)_capilares', folder_name)

                islet = islet_match.group(1) if islet_match else None
                glucose = glucose_match.group(1) if glucose_match else None
                capillary = int(capillary_match.group(1)) if capillary_match else None

                # Apply filters
                if islets and islet not in islets:
                    continue
                if glucose_levels and glucose not in glucose_levels:
                    continue
                if capillaries and capillary not in capillaries:
                    continue

                # Extract data from the H5*_data_test.dat file
                try:
                    with open(file_path, 'r') as f:
                        lines = f.readlines()
                        if lines:
                            last_line = lines[-1].strip()
                            columns = last_line.split()
                            avg_total_oxygen = float(columns[1])  # Second column
                            avg_medium_oxygen = float(columns[2])  # Third column
                        else:
                            avg_total_oxygen = None
                            avg_medium_oxygen = None
                except Exception as e:
                    print(f"Error reading file {file_path}: {e}")
                    avg_total_oxygen = None
                    avg_medium_oxygen = None

                # Correct the file_name by removing the '_data_test' suffix
                corrected_file_name = file_name.replace('_data_test', '')

                # Construct the path for the H5*_estado_cells_test.dat file
                cell_file_path = os.path.join(root, f"{corrected_file_name}_estado_cells_test.dat")
                print(f"Checking file: {cell_file_path}")
                cell_file_path = os.path.normpath(cell_file_path)  # Normalize the path
                if os.path.isfile(cell_file_path):  # Check if the file exists
                    print(f"File exists: {cell_file_path}")
                else:
                    print(f"File does not exist: {cell_file_path}")  # Debugging message for missing file
                num_cells = 0
                num_alpha = 0
                num_beta = 0
                num_delta = 0
                avg_cell_oxygen = 0
                alpha_oxygen = []
                beta_oxygen = []
                delta_oxygen = []
                functional = 0
                hypoxic = 0
                non_viable = 0

                # Initialize variables with default values
                avg_alpha_oxygen = None
                avg_beta_oxygen = None
                avg_delta_oxygen = None
                functional_prop = None
                hypoxic_prop = None
                non_viable_prop = None

                # Initialize variables for proportions by population
                functional_alpha_prop = None
                hypoxic_alpha_prop = None
                non_viable_alpha_prop = None

                functional_beta_prop = None
                hypoxic_beta_prop = None
                non_viable_beta_prop = None

                functional_delta_prop = None
                hypoxic_delta_prop = None
                non_viable_delta_prop = None

                # Initialize counters for functional, hypoxic, and non-viable cells by population
                functional_alpha = 0
                hypoxic_alpha = 0
                non_viable_alpha = 0

                functional_beta = 0
                hypoxic_beta = 0
                non_viable_beta = 0

                functional_delta = 0
                hypoxic_delta = 0
                non_viable_delta = 0

                if os.path.isfile(cell_file_path):
                    #print("success")
                    try:
                        with open(cell_file_path, 'r') as f:
                            print(cell_file_path)
                            lines = f.readlines()
                            num_cells = len(lines)
                            for line in lines:
                                #print(f"Line: {line.strip()}")
                                try:
                                    parts = line.split()
                                    if len(parts) != 4:
                                        print(f"Skipping malformed line: {line.strip()}")
                                        continue
                                    _, cell_type, oxygen, state = map(float, parts)
                                    #print(f"Parsed - Cell Type: {cell_type}, Oxygen: {oxygen}, State: {state}")
                                    oxygen = float(oxygen)
                                    state = int(state)

                                    # Count cell types
                                    if cell_type == 2:
                                        num_alpha += 1
                                        alpha_oxygen.append(oxygen)
                                        if state == 1:
                                            functional_alpha += 1
                                        elif state == 2:
                                            hypoxic_alpha += 1
                                        elif state == 3:
                                            non_viable_alpha += 1
                                    elif cell_type == 1:
                                        num_beta += 1
                                        beta_oxygen.append(oxygen)
                                        if state == 1:
                                            functional_beta += 1
                                        elif state == 2:
                                            hypoxic_beta += 1
                                        elif state == 3:
                                            non_viable_beta += 1
                                    elif cell_type == 3:
                                        num_delta += 1
                                        delta_oxygen.append(oxygen)
                                        if state == 1:
                                            functional_delta += 1
                                        elif state == 2:
                                            hypoxic_delta += 1
                                        elif state == 3:
                                            non_viable_delta += 1

                                    # Count cell states
                                    if state == 1:
                                        functional += 1
                                    elif state == 2:
                                        hypoxic += 1
                                    elif state == 3:
                                        non_viable += 1
                                except ValueError as ve:
                                    print(f"Error parsing line: {line.strip()} - {ve}")
                                    continue

                            # Calculate averages
                            avg_cell_oxygen = sum(alpha_oxygen + beta_oxygen + delta_oxygen) / num_cells if num_cells > 0 else 0
                            avg_alpha_oxygen = sum(alpha_oxygen) / len(alpha_oxygen) if alpha_oxygen else 0
                            avg_beta_oxygen = sum(beta_oxygen) / len(beta_oxygen) if beta_oxygen else 0
                            avg_delta_oxygen = sum(delta_oxygen) / len(delta_oxygen) if delta_oxygen else 0

                            # Calculate proportions
                            functional_prop = functional / num_cells if num_cells > 0 else 0
                            hypoxic_prop = hypoxic / num_cells if num_cells > 0 else 0
                            non_viable_prop = non_viable / num_cells if num_cells > 0 else 0

                            # Calculate proportions for alpha cells
                            functional_alpha_prop = functional_alpha / num_alpha if num_alpha > 0 else 0
                            hypoxic_alpha_prop = hypoxic_alpha / num_alpha if num_alpha > 0 else 0
                            non_viable_alpha_prop = non_viable_alpha / num_alpha if num_alpha > 0 else 0

                            # Calculate proportions for beta cells
                            functional_beta_prop = functional_beta / num_beta if num_beta > 0 else 0
                            hypoxic_beta_prop = hypoxic_beta / num_beta if num_beta > 0 else 0
                            non_viable_beta_prop = non_viable_beta / num_beta if num_beta > 0 else 0

                            # Calculate proportions for delta cells
                            functional_delta_prop = functional_delta / num_delta if num_delta > 0 else 0
                            hypoxic_delta_prop = hypoxic_delta / num_delta if num_delta > 0 else 0
                            non_viable_delta_prop = non_viable_delta / num_delta if num_delta > 0 else 0
                    except Exception as e:
                        print(f"Error reading file {cell_file_path}: {e}")

                # Append results
                results.append({
                    # 'File': file_name,  # Commented out to remove the File column
                    'Folder': folder_name,
                    'Islet': islet,
                    'Glucose': glucose,
                    'Capillaries': capillary,
                    'Avg_Total_Oxygen': avg_total_oxygen,
                    'Avg_Medium_Oxygen': avg_medium_oxygen,
                    'Num_Cells': num_cells,
                    'Num_Alpha': num_alpha,
                    'Num_Beta': num_beta,
                    'Num_Delta': num_delta,
                    'Avg_Cell_Oxygen': avg_cell_oxygen,
                    'Avg_Alpha_Oxygen': avg_alpha_oxygen,
                    'Avg_Beta_Oxygen': avg_beta_oxygen,
                    'Avg_Delta_Oxygen': avg_delta_oxygen,
                    'Functional_Prop': functional_prop,
                    'Hypoxic_Prop': hypoxic_prop,
                    'Non_Viable_Prop': non_viable_prop,
                    'Functional_Alpha_Prop': functional_alpha_prop,
                    'Hypoxic_Alpha_Prop': hypoxic_alpha_prop,
                    'Non_Viable_Alpha_Prop': non_viable_alpha_prop,
                    'Functional_Beta_Prop': functional_beta_prop,
                    'Hypoxic_Beta_Prop': hypoxic_beta_prop,
                    'Non_Viable_Beta_Prop': non_viable_beta_prop,
                    'Functional_Delta_Prop': functional_delta_prop,
                    'Hypoxic_Delta_Prop': hypoxic_delta_prop,
                    'Non_Viable_Delta_Prop': non_viable_delta_prop,
                })

    # Convert results to a DataFrame
    results_df = pd.DataFrame(results)
    return results_df


if __name__ == "__main__":
    import argparse

    # Set up argument parser
    parser = argparse.ArgumentParser(description="Analyze .dat files in folders based on user-defined criteria.")
    parser.add_argument("base_folder", type=str, help="Base folder containing subfolders with data files.")
    parser.add_argument("--islets", nargs="+", help="List of islet identifiers to filter (e.g., H51 H52).")
    parser.add_argument("--glucose_levels", nargs="+", help="List of glucose levels to filter (e.g., G6 G7).")
    parser.add_argument("--capillaries", nargs="+", type=int, help="List of capillary counts to filter (e.g., 5 10).")

    args = parser.parse_args()

    # Run the analysis
    results_df = analyze_files(
        base_folder=args.base_folder,
        islets=args.islets,
        glucose_levels=args.glucose_levels,
        capillaries=args.capillaries
    )

    # Save results to a CSV file
    output_file = os.path.join(args.base_folder, "analysis_results.csv")
    results_df.to_csv(output_file, index=False)
    print(f"Analysis complete. Results saved to {output_file}.")