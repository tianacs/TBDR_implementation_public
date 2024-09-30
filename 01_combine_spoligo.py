import pandas as pd
import os
import json

def combine_data_from_directory(source_folder):
    """
    Extracts data from all JSON files in the specified directory and saves it to a CSV file.

    Args:
        source_folder (str): The path to the directory containing the JSON files.

    Returns:
        pandas.DataFrame: A DataFrame containing the extracted data from all JSON files.
    """
    df = pd.DataFrame()

    # Walk through the subfolders in the source folder
    for root, dirs, files in os.walk(source_folder):
        for file in files:
            if file.endswith('results.json'):
                # Get the name of the directory the file is in
                directory_name = os.path.basename(root)
                
                # Construct full file path     
                file_path = os.path.join(source_folder, directory_name, file)
                with open(file_path, 'r') as file:
                    data = json.load(file)
                    # Process the data from the JSON file
                    df = pd.concat([df, extract_data(data, directory_name)])
    return df

def extract_data(data, name):
    """
    Extracts specific data from a JSON object and returns it as a pandas DataFrame.

    Args:
        data (dict): The JSON object containing the data.
        name (str): The name of the experiment.

    Returns:
        pandas.DataFrame: The extracted data in the form of a DataFrame.

    Raises:
        KeyError: If the required keys are not present in the JSON object.

    """
    extracted_data = []
    for sample in data['samples']:
        lineage = None
        speciesID = None
        speciesTax = None
        if sample['sample_type'] == 'test_sample':
            # Extract the lineage from the JSON data
            if sample['results']['spoligotype']['sitvit2']:
                lineage = sample['results']['spoligotype']['sitvit2'][0]['lineage']
            else:
                lineage = None
            # Extract the species identified
            if sample['results']['species_identification']['detected_species']:
                speciesTax = sample['results']['species_identification']['detected_species'][0]['taxonomic_id']
                speciesID = sample['results']['species_identification']['detected_species'][0]['scientific_name']
            extracted_data.append({
                'sample_id': sample['alias'],
                'sample_type': sample['sample_type'],
                'experiment': name,
                'sample_pass': sample['sample_pass'],  # did the sample pass quality checks?
                'sample_coverage': sample['sample_checks'][0]['check_pass'], # did the sample pass quality check for AMR
                'spoligotyping': sample['results']['spoligotype']['octal'],
                'lineage': lineage,
                'spoligotyping_qc': sample['results']['spoligotype']['spoligotype_status'], # did the sample pass quality check for spoligotyping
                'speciesTax' : speciesTax,
                'speciesID' : speciesID,
            }) 
    df = pd.DataFrame(extracted_data)
    return df

# Example usage
directory_path = r'D:\Sequencing data\Results from updated pipeline (IFIK)\Output\IeDEA_reanalyzed_in_CH_output'
df = combine_data_from_directory(directory_path)
df.to_csv('data/wf_tb_amr_v2.0.0-alpha4_spoligo.csv', index=False)
