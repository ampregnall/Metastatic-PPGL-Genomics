import os
import sys
import pandas as pd


def combine_csv_files(input_folder, rownames_file, output_file):
    # Read the row names file
    rownames_df = pd.read_csv(rownames_file)

    # Initialize an empty DataFrame to store the combined data
    combined_data = pd.DataFrame()

    # Iterate through each CSV file and concatenate its data
    for file_name in os.listdir(input_folder):
        if file_name.endswith("corrected_profile.csv"):
            file_path = os.path.join(input_folder, file_name)
            df = pd.read_csv(file_path)
            combined_data = pd.concat([combined_data, df], axis=1)

    # Combine row names with the concatenated data
    result = pd.concat([rownames_df, combined_data], axis=1)

    # Write the result to a new CSV file
    result.to_csv(output_file, index=False, sep=",")


if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python script.py input_folder rownames_file output_file")
        sys.exit(1)

    # Get command-line arguments
    input_folder = sys.argv[1]
    rownames_file = sys.argv[2]
    output_file = sys.argv[3]

    combine_csv_files(input_folder, rownames_file, output_file)
