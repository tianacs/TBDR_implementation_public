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
    ntc_control = data['workflow_checks'][0]['check_pass']
    pc_control = data['workflow_checks'][1]['check_pass']
    for sample in data['samples']:
        if sample['sample_type'] == 'test_sample':
            extracted_data.append({
                'sample_id': sample['alias'],
                #'sample_type': sample['sample_type'],
                'experiment': name,
                'ntc_control': ntc_control,
                'pc_control': pc_control,
                'sample_pass': sample['sample_pass'],  # did the sample pass quality checks?
                'sample_coverage': sample['sample_checks'][0]['check_pass'], # did the sample pass quality check for AMR
                'hsp_coverage': sample['sample_checks'][1]['check_pass'],
                'ic_coverage' : sample['sample_checks'][2]['check_pass'],
                'AMK': sample['results']['antimicrobial_resistance']['AMK']['antimicrobial_resistance_status'],
                'BDQ': sample['results']['antimicrobial_resistance']['BDQ']['antimicrobial_resistance_status'],
                'CAP': sample['results']['antimicrobial_resistance']['CAP']['antimicrobial_resistance_status'],
                'CFZ': sample['results']['antimicrobial_resistance']['CFZ']['antimicrobial_resistance_status'], 
                'DLM': sample['results']['antimicrobial_resistance']['DLM']['antimicrobial_resistance_status'], 
                'EMB': sample['results']['antimicrobial_resistance']['EMB']['antimicrobial_resistance_status'],
                'ETO': sample['results']['antimicrobial_resistance']['ETO']['antimicrobial_resistance_status'], 
                'INH': sample['results']['antimicrobial_resistance']['INH']['antimicrobial_resistance_status'], 
                'KAN': sample['results']['antimicrobial_resistance']['KAN']['antimicrobial_resistance_status'], 
                'LFX': sample['results']['antimicrobial_resistance']['LFX']['antimicrobial_resistance_status'], 
                'LZD': sample['results']['antimicrobial_resistance']['LZD']['antimicrobial_resistance_status'], 
                'MXF': sample['results']['antimicrobial_resistance']['MXF']['antimicrobial_resistance_status'], 
                'PMD': sample['results']['antimicrobial_resistance']['PMD']['antimicrobial_resistance_status'], 
                'PZA': sample['results']['antimicrobial_resistance']['PZA']['antimicrobial_resistance_status'], 
                'RIF': sample['results']['antimicrobial_resistance']['RIF']['antimicrobial_resistance_status'], 
                'STM': sample['results']['antimicrobial_resistance']['STM']['antimicrobial_resistance_status']
            }) 
    df = pd.DataFrame(extracted_data)
    return df

# Example usage
#directory_path = r'E:\Sequencing data\output'
directory_path = r'D:\Sequencing data\Results from updated pipeline (IFIK)\Output\IeDEA_reanalyzed_in_CH_output'
df = combine_data_from_directory(directory_path)
df.to_csv('data/wf_tb_amr_v2.0.0-alpha4_json.csv', index=False)
