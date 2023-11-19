import csv
import argparse


def select_column(input_file, output_file, column_name):
    with open(input_file, "r") as csv_file:
        reader = csv.DictReader(csv_file)
        headers = reader.fieldnames

        if column_name not in headers:
            print(f"Column '{column_name}' not found in the CSV file.")
            return

        selected_column = [row[column_name] for row in reader]

    with open(output_file, "w", newline="") as csv_output:
        writer = csv.writer(csv_output)
        writer.writerow([column_name])
        writer.writerows(zip(selected_column))

    print(f"Selected column '{column_name}' has been written to '{output_file}'.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select a column from a CSV file.")
    parser.add_argument("input_file", help="Path to the input CSV file")
    parser.add_argument("output_file", help="Path to the output CSV file")
    parser.add_argument("column_name", help="Name of the column to select")

    args = parser.parse_args()

    select_column(args.input_file, args.output_file, args.column_name)
