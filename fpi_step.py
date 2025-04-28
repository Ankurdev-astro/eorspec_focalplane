"""
This script calculates annulus information for a given FPI step and band type.
It reads radius and frequency data from a CSV file, performs calculations, 
and writes the results to a text file. The script loops through 20 annuli 
for each band type and handles errors gracefully. It stops when empty annuli/ freq info is hit.
"""

import pandas as pd
import numpy as np
import csv, os

def calculate_annulus(csv_file, column_name, annulus_num):
    """
    Calculates the annulus information for a given annulus number.
    
    Arguments:
    ----------
    csv_file : str
        The path to the CSV file containing the radius and frequency data.
    
    column_name : str
        The column name containing the radius data for the specified FPI step and band type.
        Format: r_step210_LFA
    
    annulus_num : int
        The annulus number for which to calculate the information.
        
    Returns:
    --------
    tuple or str
        Returns a tuple containing the annulus information if successful.
        The tuple contains the following values:
        r_min, r_max, freq_min, freq_max, freq_centre, freq_delta
        
        Returns an error message if an issue occurs during the calculation.        
    
    """
    
    row_num = annulus_num - 1 # For zero indexing
    if row_num < 0:
        raise ValueError('Minimum Value for Annulus is 0')
    try:
        # Reading Data
        df = pd.read_csv(csv_file)

        column_start = f"{column_name}_start"
        column_end = f"{column_name}_end"
        # Data Cleaning: Drop NaN rows for the specified column
        df_cleaned = df.dropna(subset=[column_start])

        # Data Extraction: Extract radius and frequency values
        radius_start = df_cleaned[column_start].values
        radius_end = df_cleaned[column_end].values
        frequency_min = df_cleaned['frequency_start[GHz]'].values
        frequency_max = df_cleaned['frequency_end[GHz]'].values

        # Annulus Calculation: Determine r_min, r_max, freq_min, and freq_max
        r_min = round(radius_start[row_num], 18)
        r_max = round(radius_end[row_num], 18)

        freq_min = frequency_min[row_num]
        freq_max = frequency_max[row_num]

        freq_delta = freq_max - freq_min
        freq_centre = (freq_max + freq_min)/2.0

        return r_min, r_max, round(freq_min,2), round(freq_max,2), \
                round(freq_centre,2), round(freq_delta,2)
    
    except FileNotFoundError:
        return "Error: CSV file not found."
    except KeyError:
        return f"Error: Column {column_name} not found."
    except IndexError:
        return f"Error: No frequency data for Annulus A{annulus_num}."
    
    
def annulus_FPIstep(step, csv_file = 'annulus_radii.csv'):
    """
    Generates a text file containing annulus information for both LFA and HFA bands for a given FPI step.
    
    Parameters:
    -----------
    step : str
        The FPI step in the format "stepXXX" (e.g., "step210").
    
    csv_file : str, optional
        The path to the CSV file containing the radius and frequency data. 
        Default is 'annulus_radii.csv'.
    
    Returns: None
    --------
    File name str.
        Writes the annulus information to a text file named `annulus_results_{step}.txt`.
        Each line in the text file contains tab-separated values in the following format:
        band_type, annulus_num, r_min, r_max, freq_min, freq_max, freq_centre, freq_delta
    
    Notes:
    ------
    - The function loops through 20 annuli for each band type ("LFA" and "HFA").
    - The function will print an error message if an issue occurs during the calculation for any annulus.
    
    Examples:
    ---------
    >>> annulus_FPIstep("step210")
    """
    
    # Define LFA and HFA filter ranges
    LFA_min_freq = 209
    LFA_max_freq = 316
    
    HFA_min_freq = 315
    HFA_max_freq = 422
    
    # Create the directory for storing annuli data if it doesn't exist
    annuli_dir = './fpi_data/annuli_data/'
    if not os.path.exists(annuli_dir):
        os.makedirs(annuli_dir)
        print(f"Created directory: {annuli_dir}")
    
    # Initialize the text file
    f_write = f"./fpi_data/annuli_data/annulus_results_{step}.txt"

    with open(f_write, "w") as f:
        f.write("wtype\tannulus_num\tfreq_channel\tr_min\tr_max \tfreq_min\tfreq_max\tfreq_centre\tfreq_delta\n")

    # Loop over the band types (LFA and HFA)
    for wtype in ["LFA", "HFA"]:
        column_name = f"r_{step}_{wtype}"
        # Loop over the annulus numbers (up to 20 as specified)
        for annulus_num in range(1, 21):  # Adjust the range according to your needs
            result = calculate_annulus(csv_file, column_name, annulus_num)
            if isinstance(result, tuple):
                r_min, r_max, freq_start, freq_end, freq_centre, freq_delta = result
                freq_channel = int(np.floor(freq_centre))
                
                # Additional conditions based on band type and Freq
                # if (wtype == "LFA" and LFA_min_freq <= freq_start and freq_end <= LFA_max_freq) or \
                #    (wtype == "HFA" and HFA_min_freq <= freq_start and freq_end <= HFA_max_freq):
                if (wtype == "LFA" and freq_end >= LFA_min_freq and freq_start <= LFA_max_freq) or \
                    (wtype == "HFA" and freq_end >= HFA_min_freq and freq_start <= HFA_max_freq):
                    
                    # Write to the text file
                    with open(f_write, "a") as f:
                        f.write(
                                f"{wtype}\t{annulus_num}\t{freq_channel}\t{r_min}\t{r_max}\t"
                                f"{freq_start}\t{freq_end}\t{freq_centre}\t{freq_delta}\n"
                                )
                else:
                    print(f"Skipping annulus {annulus_num} for {wtype} due to frequency range.")
                    print(f'Wtype: {wtype}, Freq min: {freq_start}, Freq max: {freq_end}', '\n')


            elif isinstance(result, str):
                print(f"An error occurred for HFA and annulus {annulus_num}: {result}")
                break
    return None


def process_steps(csv_file='annulus_radii.csv'):
    """
    Processes the annulus_radii CSV file and returns a list of FPI steps.
    
    Parameters:
    -----------
    csv_file : str, optional
        The path to the CSV file containing the radius and frequency data. 
        Default is 'annulus_radii.csv'.
        
    Returns:
    --------
    list
        A list of FPI steps extracted from the CSV file.
    """
    
    fpi_steps = []

    with open(csv_file, mode='r') as file:
        # reading the CSV file
        csvFile = csv.reader(file)

        # displaying the contents of the CSV file
        for lines in csvFile:
            column_names = lines
            break  # Stop after first line

    # Filter out the relevant column names and populate the fpi_steps dictionary
    for name in column_names:
        if "_LFA_start" in name:
            step = name.split('_')[1]  # Assuming the format is always r_stepXXX_LFA_start
            fpi_steps.append(step)
    
    return fpi_steps


def getall_freq_chl():
    """
    Returns a list of all unique frequency channels from the annulus results files.
       
    """
    
    # Get the list of fpi steps
    fpi_steps = process_steps()
    
    # Base directory for storing HDF5 files
    annuli_dir = './fpi_data/annuli_data/'
    freq_channel_list = []
    for i, step in enumerate(fpi_steps):

        fpistep_infotxt = os.path.join(annuli_dir, f"annulus_results_{step}.txt")
        # Read the annulus results file
        df = pd.read_csv(fpistep_infotxt, delimiter='\t')
        
        # Get the unique freq_channel values
        freq_channels = df['freq_channel'].unique()
        
        # Append the freq_channel values to the list if they don't exist already
        for freq_chl in freq_channels:
            if freq_chl not in freq_channel_list:
                freq_channel_list.append(freq_chl)
        
    # print(f"Freq channel list: {freq_channel_list}")
    return freq_channel_list

if __name__ == "__main__":
    fpi_steps = process_steps()
    print(fpi_steps)

    print("\n","="*50,"\n")
    for step in fpi_steps:
        print("\n","-"*50,"\n")
        print(f"Processing FPI step {step}...")
        annulus_FPIstep(step)