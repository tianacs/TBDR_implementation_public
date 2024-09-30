import glob
import os
import pandas as pd

def combine_csv_files(source_folder):
    """
    Combine multiple CSV files into a single DataFrame.

    Args:
        folder_path (str): The path to the root folder containing the CSV files.

    Returns:
        pandas.DataFrame: The combined DataFrame containing data from all CSV files.
    """
    # Generate the list of all .csv files in the root folder and its subfolders
    csv_files = glob.glob(os.path.join(source_folder, '**', '*.csv'), recursive=True)
    
    # Read each CSV file into a DataFrame and store in a list
    dataframes = []
    for file in csv_files:
        df = pd.read_csv(file)
        # Add folder name as a column in the DataFrame
        folder_name = os.path.basename(os.path.dirname(file))
        df['experiment'] = folder_name
        dataframes.append(df)

    # Combine all DataFrames into a single DataFrame
    combined_df = pd.concat(dataframes, ignore_index=True)

    # Print the combined DataFrame (optional)
    return combined_df

# Example usage
#source_folder = r'E:\Sequencing data\output'
source_folder = r'D:\Sequencing data\Results from updated pipeline (IFIK)\Output\IeDEA_reanalyzed_in_CH_output'
df = combine_csv_files(source_folder)

# Save the combined DataFrame to a CSV file
df.to_csv('data/wf_tb_amr_v2.0.0-alpha4_csv.csv', index=False)
