### PACKAGES ###

import os
import re
import time
import scipy
import pickle
import shutil
import logging
import sklearn
import warnings
import statistics
import subprocess
import numpy as np
import Bio.SeqUtils
import pandas as pd
import seaborn as sns
import plotly.express as px
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from itertools import combinations
from sklearn.cluster import KMeans
from scipy.spatial import distance
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
from Bio.PDB.PDBParser import PDBParser
from tempfile import NamedTemporaryFile
from scipy.spatial.distance import cdist
from scipy.spatial.distance import pdist
from scipy.spatial.distance import euclidean
from Bio.PDB.Superimposer import Superimposer
from sklearn.preprocessing import MinMaxScaler
from matplotlib.collections import PolyCollection
from scipy.spatial.distance import pdist, squareform
from subprocess import CalledProcessError, run, PIPE
from Bio.PDB import PDBParser, MMCIFParser, Selection
from sklearn.model_selection import StratifiedShuffleSplit
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from itertools import combinations, combinations_with_replacement, permutations
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

warnings.filterwarnings('ignore') # ignores warnings
logging.disable(logging.CRITICAL) # disables logging
log = logging.getLogger("LIGYSIS") # starts ligysis logger

from prointvar.fetchers import download_structure_from_pdbe
from prointvar.pdbx import PDBXreader, PDBXwriter
from prointvar.config import config as cfg

### FUNCTIONS ###

def save_to_pickle(variable, file_path):
    with open(file_path, 'wb') as file:
        pickle.dump(variable, file)

def read_from_pickle(file_path):
    with open(file_path, 'rb') as file:
        return pickle.load(file)

def download_and_move_files(pdb_ids, asymmetric_dir, bio = False):
    """
    Downloads CIF of a series of PDB IDs and moves
    them to a given directory.
    """
    cifs = []
    for pdb_id in pdb_ids:
        if bio:
            cif_in = os.path.join(cfg.db_root, cfg.db_pdbx, "{}_bio.cif".format(pdb_id))
            cif_out = os.path.join(asymmetric_dir, "{}_bio.cif".format(pdb_id))
        else:
            cif_in = os.path.join(cfg.db_root, cfg.db_pdbx, "{}.cif".format(pdb_id))
            cif_out = os.path.join(asymmetric_dir, "{}.cif".format(pdb_id))
        if os.path.isfile(cif_out):
            log.debug("{} already exists!".format(cif_out))
        else:
            download_structure_from_pdbe(pdb_id, bio = bio)
            shutil.move(cif_in, cif_out)
        cifs.append(cif_out)
    return cifs

def fix_col_names(df):
	"""
	Originally written to fix P2Rank output
	dataframe columns, which contained white
	spaces.
	"""
	cols = df.columns.tolist()
	new_cols = [col.strip() for col in cols]
	df.columns = new_cols
	return df

def reformat_sites_df(df):
    """
    Reformats P2Rank predicted pockets dataframe.
    """
    # Copy the DataFrame to avoid modifying the original data
    df_transformed = df.copy()
    
    # Step 1: Rename 'name' column
    df_transformed['name'] = df_transformed['name'].str.replace('pocket', '').astype(int)
    
    # Step 2: Convert 'residue_ids' to list of integers
    df_transformed['residue_ids'] = df_transformed['residue_ids'].apply(
        lambda x: [id_.split("_")[-1] for id_ in x.split()]
    )

    df_transformed['centre'] = list(zip(df_transformed['center_x'], df_transformed['center_y'], df_transformed['center_z']))
    
    df_transformed.drop(columns=["center_x", "center_y", "center_z"], inplace = True)
    
    df_transformed['n_aas'] = df_transformed['residue_ids'].apply(lambda x: len(x)) 
    
    # Step 3: Rename columns
    column_renames = {
        'name': 'ID',
        'rank': 'RANK',
        'score': 'score',
        'probability': 'prob',
        'sas_points': 'n_sas_points',
        'surf_atoms': 'n_surf_atoms',
        'residue_ids': 'aas'
    }
    df_transformed.rename(columns=column_renames, inplace=True)
    
    # Add a new column for site size

    df_transformed.n_sas_points = df_transformed.n_sas_points.astype(int)
    df_transformed.n_surf_atoms = df_transformed.n_surf_atoms.astype(int)
    
    return df_transformed

def read_files_to_dataframe(directory):
    """
    This function parses all the freeSASA output
    files (whole protein) found in a given directory
    and returns a dataframe.
    """
    # Initialize a list to collect data
    data = []

    # List all files in the given directory
    for filename in os.listdir(directory):
        if filename.endswith('.out'):  # adjust to '.out' as per your file extension
            file_path = os.path.join(directory, filename)
            with open(file_path, 'r') as file:
                contents = file.read()
            
            # Parse the file contents
            parsed_data = parse_freesasa_results(contents)
            
            # Create a temporary dictionary to hold results for this file
            temp_data = {'rep_chain': filename.split(".")[0]}
            results_data = parsed_data.get('RESULTS', {})

            # Add results to the temporary dictionary
            for key in results_data:
                temp_data[key] = results_data[key]

            # Append the temporary dictionary to the data list
            data.append(temp_data)

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(data)
    return df

def calculate_rmsd(file1, file2):
    # Initialize parser
    parser = PDBParser()

    # Read structures from files
    structure1 = parser.get_structure('struc1', file1)
    structure2 = parser.get_structure('struc2', file2)

    # Extract first model
    model1 = structure1[0]
    model2 = structure2[0]

    # Select atoms: Here we are using CA atoms (can change based on your molecule)
    atoms1 = [atom for atom in model1.get_atoms() if atom.name == 'CA']
    atoms2 = [atom for atom in model2.get_atoms() if atom.name == 'CA']

    # Ensure same number of atoms
    assert len(atoms1) == len(atoms2), "ERROR: Structures do not have the same number of CA atoms."

    # Superimpose the atom sets using Superimposer
    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    sup.apply(atoms2)  # This applies the rotation/translation to the second set of atoms

    # Return the RMSD
    return sup.rms

def extract_volume_data(directory):
    """
    This function parses all the Protein Volume output
    files in a directory and returns a dataframe.
    """
    records = []  # List to hold all records

    for filename in os.listdir(directory):
        if filename.endswith(".protvol.out"):  # Adjusted for the actual file extension
            filepath = os.path.join(directory, filename)
            with open(filepath, 'r') as file:
                content = file.readlines()
                
                protein_path = 'Unknown'
                for line in content:
                    if line.startswith('Protein:'):
                        protein_path = line.split('Protein:')[1].strip()
                        protein_name = os.path.basename(protein_path)
                
                # Directly look for the data line
                for line in content:
                    if line.strip() and not line.startswith(('Could not', 'Using', 'PV Version', 'Date', 'Protein', 'Volume probe', 'Surface probe', 'Reading', 'Number of')):
                        data_line = line.strip()
                        break
                
                # Parse the data line
                try:
                    _, total_volume, void_volume, vdw_volume, packing_density, time_taken = data_line.split(maxsplit=5)
                    records.append({
                        "rep_chain": protein_name.split(".")[0],
                        "total_volume": float(total_volume),
                        "void_volume": float(void_volume),
                        "vdw_volume": float(vdw_volume),
                        "packing_density": float(packing_density),
                        "time_taken": int(time_taken)
                    })
                except ValueError:
                    # Log the error or handle files that do not match the expected format
                    print(f"Error parsing file: {filepath}")

    df = pd.DataFrame(records)
    return df

def extract_residue_number(residue_info):
    """
    Extract integer part of the residue number (ignoring any non-digit characters like '84F')
    """
    residue_number = ''.join(filter(str.isdigit, residue_info))
    return int(residue_number) if residue_number else None

def parse_breaks_file(filepath):
    """
    Function to parse the '.breaks' file
    """
    with open(filepath, 'r') as file:
        #print(file.readlines())
        breaks_data = file.readlines()
    breaks_data = [eval(break_line.strip()) for break_line in breaks_data]
    return breaks_data

def parse_break_residues_file(filepath):
    """
    Function to parse the '.break_residues' file
    """
    with open(filepath, 'r') as file:
        break_residues_data = file.readlines()
    break_residues_data = [tuple(res.strip().split(',')) for res in break_residues_data]
    return break_residues_data

def find_center_and_dimensions(coords, buffer = 5):
    """
    This function determines an inclusion region to calculate
    POVME and calculate a pocket volume. Takes as input the 3D
    coordinates of the atom forming the site.
    """
    xs, ys, zs = zip(*coords)
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)
    
    # Calculate center
    center_x = round((max_x + min_x) / 2, 2)
    center_y = round((max_y + min_y) / 2, 2)
    center_z = round((max_z + min_z) / 2, 2)
    
    # Calculate dimensions
    width_x = round(max_x - min_x + buffer, 2) 
    width_y = round(max_y - min_y + buffer, 2) 
    width_z = round(max_z - min_z + buffer, 2) 
    
    return (center_x, center_y, center_z, width_x, width_y, width_z)

def run_freesasa(probe_radius, n_threads, input_file, output_file, error_file, format = 'seq', n_slices = 20, override = False):
    """
    Run the FreeSASA command with specified parameters.

    Args:
    probe_radius (float): The probe radius to use.
    n_threads (int): The number of threads to use.
    input_file (str): The path to the input PDB file.
    output_file (str): The path to the output file.
    error_file (str): The path to the error file.
    """

    if os.path.isfile(output_file) and override == False:
        #print(f'{output_file} already exists!')
        return
    
    command = [
        'freesasa',
        '--lee-richards',
        '--probe-radius', str(probe_radius),
        '--n-threads=' + str(n_threads),
        '--format=' + format,
        '-n', str(n_slices),
        '<', input_file,
        '--output=' + output_file,
        '--error-file=' + error_file
    ]
    
    # Execute the command
    result = subprocess.run(' '.join(command), shell=True, capture_output=True, text=True)
    #print(' '.join(command))
    
    # Check if the command was successful
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
    else:
        pass

def parse_freesasa_output(filename):
    """
    This function parses freeSASA results for
    a protein ran in 'seq' mode where you get
    the individual SASA per residue. It also
    normalises it using the method presented
    by Tien et al., 2013.
	"""
    sasa_dict = {}
    rsa_dict = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('SEQ'):
                parts = line.split()
                chain_id = parts[1]
                res_num = parts[2]
                res_name = parts[3]
                sasa = float(parts[5])
                #key = (chain_id, res_num, res_name)
                sasa_dict[res_num] = sasa
                try:
                    RSA = round(100*(sasa/ASA_Wilke[res_name]), 2)
                except:
                    RSA = round(100*(sasa/avg_max_sasa), 2)
                rsa_dict[res_num] = RSA
    return sasa_dict, rsa_dict

def parse_freesasa_results(contents):
    """
    This function parses freeSASA results for
    a whole protein.
    """
    lines = contents.strip().split('\n')
    parsed_data = {}
    current_block = None

    for line in lines:
        line = line.strip()
        if not line or line.startswith('##') or line.endswith('##'):
            continue
        if line.isupper():
            current_block = line.replace(' (A^2)', '')
            parsed_data[current_block] = {}
        else:
            key, value = line.split(':')
            key = key.strip()
            value = value.strip()
            try:
                if '.' in value:
                    value = float(value)
                else:
                    value = int(value)
            except ValueError:
                pass
            parsed_data[current_block][key] = value

    return parsed_data

def run_povme(command, site_id, retries=3, delay=2):
    """
    Runs POVME given a command, an identifier for a binding site,
    the number of retries and delay between them.
    """
    for attempt in range(retries):
        try:
            result = run(command, check=True, stdout=PIPE, stderr=PIPE)
            return result
        except CalledProcessError as e:
            print(f"Attempt {attempt+1} failed for {site_id}")
            time.sleep(delay)  # wait before retrying
    raise Exception("All attempts to run POVME failed.")

def calculate_total_sasa(df, sasa_dict, column = 'aas'):
    """
    Calculate the total SASA for each site in a pandas DataFrame.

    Parameters:
    - df: pandas DataFrame with columns ['origin', 'rep_chain', 'ID', 'aas']
    - sasa_dict: Dictionary with rep_chain as keys and another dict as value that maps aa positions to SASA values.

    Returns:
    - A new DataFrame with an additional 'total_sasa' column representing the sum of SASA values for each site.
    """

    # Function to sum the SASA for each amino acid in the list for a given rep_chain
    def sum_sasa(row):
        chain_sasa = sasa_dict.get(row['rep_chain'], {})
        if isinstance(row[column], list):
            return sum(chain_sasa.get(str(aa), 0) for aa in row[column])
        elif np.isnan(row[column]):
            return np.nan
        #else:
            
    
    # Apply the function to each row in the DataFrame
    df['SASA'] = df.apply(sum_sasa, axis=1)
    
    return df

def read_files_to_dataframe(directory):
    # Initialize a list to collect data
    data = []

    # List all files in the given directory
    for filename in os.listdir(directory):
        if filename.endswith('.out'):  # adjust to '.out' as per your file extension
            file_path = os.path.join(directory, filename)
            with open(file_path, 'r') as file:
                contents = file.read()
            
            # Parse the file contents
            parsed_data = parse_freesasa_results(contents)
            
            # Create a temporary dictionary to hold results for this file
            temp_data = {'rep_chain': filename.split(".")[0]}
            results_data = parsed_data.get('RESULTS', {})

            # Add results to the temporary dictionary
            for key in results_data:
                temp_data[key] = results_data[key]

            # Append the temporary dictionary to the data list
            data.append(temp_data)

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(data)
    return df

def extract_volume_data(directory):
    records = []  # List to hold all records

    for filename in os.listdir(directory):
        if filename.endswith(".protvol.out"):  # Adjusted for the actual file extension
            filepath = os.path.join(directory, filename)
            with open(filepath, 'r') as file:
                content = file.readlines()
                
                protein_path = 'Unknown'
                for line in content:
                    if line.startswith('Protein:'):
                        protein_path = line.split('Protein:')[1].strip()
                        protein_name = os.path.basename(protein_path)
                
                # Directly look for the data line
                for line in content:
                    if line.strip() and not line.startswith(('Could not', 'Using', 'PV Version', 'Date', 'Protein', 'Volume probe', 'Surface probe', 'Reading', 'Number of')):
                        data_line = line.strip()
                        break
                
                # Parse the data line
                try:
                    _, total_volume, void_volume, vdw_volume, packing_density, time_taken = data_line.split(maxsplit=5)
                    records.append({
                        "rep_chain": protein_name.split(".")[0],
                        "total_volume": float(total_volume),
                        "void_volume": float(void_volume),
                        "vdw_volume": float(vdw_volume),
                        "packing_density": float(packing_density),
                        "time_taken": int(time_taken)
                    })
                except ValueError:
                    # Log the error or handle files that do not match the expected format
                    print(f"Error parsing file: {filepath}")

    df = pd.DataFrame(records)
    return df

# Function to safely extract the integer part from a string that represents a residue number, ignoring non-digit characters
def extract_residue_number(residue_info):
    # Extract integer part of the residue number (ignoring any non-digit characters like '84F')
    residue_number = ''.join(filter(str.isdigit, residue_info))
    return int(residue_number) if residue_number else None

# Function to parse the '.breaks' file
def parse_breaks_file(filepath):
    with open(filepath, 'r') as file:
        #print(file.readlines())
        breaks_data = file.readlines()
    breaks_data = [eval(break_line.strip()) for break_line in breaks_data]
    return breaks_data

# Function to parse the '.break_residues' file
def parse_break_residues_file(filepath):
    with open(filepath, 'r') as file:
        break_residues_data = file.readlines()
    break_residues_data = [tuple(res.strip().split(',')) for res in break_residues_data]
    return break_residues_data

def extract_volumes_from_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()  # Read all lines from the file

    # Assert that there is exactly one line in the file
    assert len(lines) == 1, "File must contain exactly one line"

    # Split the line into parts and convert them to the appropriate types
    line = lines[0].strip()  # Remove any trailing whitespace
    parts = line.split()     # Split the line into parts by whitespace

    # Convert the first part to integer and the second part to float
    first_value = int(parts[0])
    second_value = float(parts[1])

    return first_value, second_value

def process_vnegnn_df(df):
    # Combine x, y, z into the "centre" column with tuples rounded to 3 significant figures
    df['centre'] = list(zip(df['x'].round(3), df['y'].round(3), df['z'].round(3)))
    
    # Rename 'rank' column to 'score'
    df.rename(columns={'rank': 'score'}, inplace=True)
    
    # Add 'ID' column which goes from 1 to N (number of rows)
    df['ID'] = range(1, len(df) + 1)
    
    # Add 'RANK' column that ranks rows by 'score' in descending order
    df['RANK'] = df['score'].rank(ascending=False, method='min').astype(int)
    
    # Drop the original x, y, z columns
    df.drop(columns=['x', 'y', 'z'], inplace=True)
    
    return df

def calculate_distance(point, centroid):
    return np.sqrt((point[0] - centroid[0]) ** 2 + (point[1] - centroid[1]) ** 2 + (point[2] - centroid[2]) ** 2)

def custom_sort(s):
    # Using regular expression to find the first group of digits at the start of the string
    # and any optional alphabetical characters that follow
    match = re.match(r"(\d+)([a-zA-Z]*)", s)
    if match:
        return int(match.group(1)), match.group(2)
    else:
        return 0, s  # Default case if no digits are found

def extract_coordinates(file_path):
    """
    Extracts X, Y, Z coordinates from a file that follows a specific format.
    
    The file format expected is:
    Line 1: a single digit (the point's ID)
    Line 2: the word "point"
    Line 3: "PS", followed by three floating-point numbers (the X, Y, Z coordinates)
    
    Parameters:
    - file_path: Path to the file to be read.
    
    Returns:
    - A tuple of floats: (X, Y, Z) coordinates.
    
    Raises:
    - AssertionError if the file does not follow the expected format.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
        # Assert there are exactly 3 lines
        assert len(lines) == 3, "File does not contain exactly three lines."
        
        # Assert the second line is exactly "point\n"
        assert lines[1].strip() == "point", "Second line is not 'point'."
        
        # Check and extract coordinates from the third line
        parts = lines[2].strip().split()
        
        # Assert the line starts with "PS" and is followed by three numbers
        assert parts[0] == "PS" and len(parts) == 4, "Third line format is incorrect."
        
        # Try converting the coordinates to floats
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            raise AssertionError("Coordinates are not valid floating-point numbers.")
        
        return x, y, z

def find_residues_near_points(structure, points, distance_threshold):
    """
    Finds residues within a given distance of specified 3D coordinates.

    Args:
    - cif_path (str): Path to the CIF file containing the structure.
    - points (list of tuples): List of 3D coordinates (x, y, z).
    - distance_threshold (float): Distance threshold for identifying nearby residues.

    Returns:
    - set: Set of residue numbers within the given distance of any of the points.
    """
    
    # Convert points to a NumPy array for efficient distance computation
    points = np.array(points)
    
    # Set to store unique residue numbers
    nearby_residues = set()
    
    # Iterate over all atoms in the structure
    for atom in Selection.unfold_entities(structure, 'A'):  # 'A' for Atoms
        atom_coord = np.array(atom.get_coord())
        
        # Calculate distances from the current atom to all points
        distances = np.linalg.norm(points - atom_coord, axis=1)
        
        # Check if any distance is within the threshold
        # if np.any(distances < distance_threshold):
        #     nearby_residues.add(atom.get_parent().get_id()[1])  # Add residue number

        if np.any(distances < distance_threshold):
            residue = atom.get_parent()
            resseq, icode = residue.get_id()[1], residue.get_id()[2]
            if icode != " ":
                print(residue.get_id()[0], residue.get_id()[1], residue.get_id()[2])
            res_id = f"{resseq}{icode.strip()}"  # Concatenate, stripping any blank spaces from icode
            nearby_residues.add(res_id)
    
    return sorted(list(nearby_residues))

def extract_multiple_coordinates(file_path):
    """
    Extracts X, Y, Z coordinates from a file that may contain multiple sets of coordinates.
    
    Each coordinate set follows a format:
    - A line with the word "point"
    - Followed by one or more lines starting with "PS", followed by three floating-point numbers (the X, Y, Z coordinates)
    
    Parameters:
    - file_path: Path to the file to be read.
    
    Returns:
    - A list of tuples, each a set of (X, Y, Z) coordinates.
    
    Raises:
    - AssertionError if any coordinate set does not follow the expected format.
    """
    coordinates = []  # Initialize an empty list to store the coordinates
    with open(file_path, 'r') as file:
        for line in file:
            # Strip newline and whitespace from the line for easier processing
            line = line.strip()
            if line.startswith("PS"):
                parts = line.split()
                # Assert the line is correctly formatted
                assert parts[0] == "PS" and len(parts) == 4, "Line format is incorrect: " + line
                try:
                    # Convert the coordinates to floats and add to the list
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    coordinates.append((x, y, z))
                except ValueError:
                    raise AssertionError("Coordinates are not valid floating-point numbers.")
    return coordinates

def apply_rotation(row):
    # Extract the rotation matrix using rep_chain
    matrix = rot_matrices[row['rep_chain']]
    # Extract the centre tuple and convert to numpy array
    try:
        centre_vector = np.array(row['centre'])
    except:
        centre_vector = np.array(row['centre_deep2'])
    translated_centre = centre_vector - orig_centroids[row['rep_chain']]
    # Perform matrix multiplication
    transformed_centre = np.dot(matrix, translated_centre)
    round_centre = tuple([round(el, 3) for el in transformed_centre])
    return round_centre

def calculate_distance(row):
    # Convert tuples to numpy arrays for vector operations
    centre_mat_vector = np.array(row['centre_mat'])
    try:
        centre_trans_vector = np.array(row['centre_trans'])
    except:
        centre_trans_vector = np.array(row['centre_deep_trans'])
    # Calculate Euclidean distance using numpy.linalg.norm
    distance = np.linalg.norm(centre_mat_vector - centre_trans_vector)
    return distance

def calculate_distance2(row):
    # Convert tuples to numpy arrays for vector operations
    centre_mat_vector = np.array(row['centre_mat'])
    try:
        centre_trans_vector = np.array(row['centre_trans'])
    except:
        centre_trans_vector = np.array(row['centre_deep_trans'])
    # Calculate Euclidean distance using numpy.linalg.norm
    distance = np.linalg.norm(centre_mat_vector - centre_trans_vector)
    return distance

def run_clean_pdb(pdb_path):
    """
    runs pdb_clean.py to prepare files for arpeggio
    """
    args = [
        clean_pdb_python_bin, clean_pdb_bin, pdb_path
    ]
    cmd = " ".join(args)
    exit_code = os.system(cmd)
    if exit_code != 0:
        print(f'ERROR with {cmd}')
    return exit_code

def parse_concavity_residue_scores(file_path):
    scores_dict = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            # Skip lines that are comments or empty
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split()
            if len(parts) == 3:
                number = int(parts[0])
                amino_acid = parts[1]
                score = float(parts[2])
                scores_dict[(number, amino_acid)] = score
    return scores_dict

def calculate_concavity_centroids(coordinates_dict):
    """
    Calculate the centroid for each set of coordinates in the input dictionary.
    
    Parameters:
    coordinates_dict (dict): Dictionary with label_seq_id_full as keys and numpy arrays of coordinates as values.
    
    Returns:
    dict: Dictionary with label_seq_id_full as keys and centroid coordinates as values.
    """
    centroids_dict = {}
    for seq_id, coordinates in coordinates_dict.items():
        centroid = np.mean(coordinates, axis=0)
        centroids_dict[seq_id] = centroid
    return centroids_dict

def create_concavity_coordinates_dict(df):
    """
    Create a dictionary where the key is label_seq_id_full and the value is a numpy array of 
    the coordinates given by Cartn_x, Cartn_y, and Cartn_z.
    
    Parameters:
    df (pd.DataFrame): Input dataframe with columns ['label_seq_id_full', 'Cartn_x', 'Cartn_y', 'Cartn_z'].
    
    Returns:
    dict: Dictionary with label_seq_id_full as keys and numpy arrays of coordinates as values.
    """
    result_dict = {}
    for seq_id, group in df.groupby('label_seq_id_full'):
        coordinates = group[['Cartn_x', 'Cartn_y', 'Cartn_z']].to_numpy()
        result_dict[seq_id] = coordinates
    return result_dict

def calculate_concavity_pocket_score(df):
    """
    Calculate the sum of squares of B_iso_or_equiv for each label_seq_id_full.
    
    Parameters:
    df (pd.DataFrame): Input dataframe with columns ['label_seq_id_full', 'B_iso_or_equiv'].
    
    Returns:
    dict: Dictionary with label_seq_id_full as keys and sum of squares of B_iso_or_equiv as values.
    """
    sum_of_squares_dict = {}
    for seq_id, group in df.groupby('label_seq_id_full'):
        sum_of_squares = np.sum(group['B_iso_or_equiv']**2)
        sum_of_squares_dict[seq_id] = sum_of_squares
    return sum_of_squares_dict

def get_concavity_pocket_ress(coords_dict, clean_df, threshold=6):
    """
    Find unique label_seq_id_full values within a given Euclidean distance threshold.
    
    Parameters:
    coords_dict (dict): Dictionary where keys are identifiers and values are np.arrays of coordinates.
    clean_df (pd.DataFrame): DataFrame containing 'Cartn_x', 'Cartn_y', 'Cartn_z', and 'label_seq_id_full'.
    threshold (float): Euclidean distance threshold.
    
    Returns:
    dict: Dictionary where keys are the same as coords_dict and values are lists of unique label_seq_id_full values.
    """
    pocket_ress = {}

    for key, arr in coords_dict.items():
        distances = cdist(arr, clean_df[['Cartn_x', 'Cartn_y', 'Cartn_z']], metric='euclidean')
        close_points = np.where(distances <= threshold)
        unique_label_seq_ids = clean_df.iloc[close_points[1]]['label_seq_id_full'].unique()
        pocket_ress[key] = sorted(unique_label_seq_ids.tolist(), key=custom_sort)
    
    return pocket_ress

def get_PRANKED_df(pranked_dir, change_id = True):
    dfs = []

    rescored_files = sorted(os.listdir(pranked_dir))
    
    for f in rescored_files:
        if f.endswith(".csv"):
            df = pd.read_csv(os.path.join(pranked_dir, f))
            df["rep_chain"] = f.split(".")[0]
            dfs.append(df)
        else:
            print(f'Error with {f}')
    pranked_df = pd.concat(dfs).reset_index(drop = True)
    
    pranked_df = pranked_df[["rep_chain", "name", "score", "rank", "old_rank", "change", "   "]].copy()

    pranked_df.rename(columns = {
        "name": "ID", "rank": "RANK", "old_rank": "OLD_RANK", "change": "CHANGE"
    }, inplace = True)

    if change_id:
        IDs = [int(el.split(".")[1]) for el in pranked_df.ID.tolist()]
        pranked_df["ID"] = IDs
    
    return pranked_df


def calculate_metrics_refined(df, origin):
    # Filter rows by origin
    filtered_df = df[df['origin'] == origin].dropna(subset=['centre', 'up_aas'])  # Drop rows where 'centre' or 'up_aas' might be NaN
    
    # Prepare a list to collect results
    results = []

    # Group by rep_chain
    for rep_chain, group in filtered_df.groupby('rep_chain'):
        # Extract data for each record
        records = [{
            'ID': row['ID'],
            'centre': np.array(row['centre']),
            'up_aas': set(row['up_aas'])
        } for index, row in group.iterrows()]
        
        # Process each unique pair of records within the group
        for rec1, rec2 in permutations(records, 2):

            # Distance calculation
            distance = np.linalg.norm(rec1['centre'] - rec2['centre'])
            
            # Set similarity metrics calculation
            intersection = rec1['up_aas'].intersection(rec2['up_aas'])
            union = rec1['up_aas'].union(rec2['up_aas'])
            jaccard_index = round(len(intersection) / len(union), 4) if union else 0
            user_relative_intersection = round(len(intersection) / min(len(rec1['up_aas']), len(rec2['up_aas'])), 4) if min(len(rec1['up_aas']), len(rec2['up_aas'])) > 0 else 0
            averaged_relative_intersection = round((len(intersection) / len(rec1['up_aas']) + len(intersection) / len(rec2['up_aas'])) / 2, 2) if len(rec1['up_aas']) and len(rec2['up_aas']) > 0 else 0
            
            # Append results
            results.append({
                'origin': origin,
                'rep_chain': rep_chain,
                'ID1': rec1['ID'],
                'ID2': rec2['ID'],
                'DIST': distance,
                'intersection': len(intersection),
                'jaccard_index': jaccard_index,
                'JSU_rel_intersection': user_relative_intersection,
                'avg_rel_intersection': averaged_relative_intersection
            })
            
    # Create a results dataframe
    results_df = pd.DataFrame(results)
    
    return results_df

def get_sites_per_prot(combined_df):
    a = combined_df.groupby('origin')['rep_chain'].value_counts()
    df = pd.DataFrame(list(a.items()), columns=['key', 'sites_per_prot'])
    df[['origin', 'rep_chain']] = pd.DataFrame(df['key'].tolist(), index=df.index)
    df.drop('key', axis=1, inplace=True)
    df = df[['origin', 'rep_chain', 'sites_per_prot']]  # Reorder columns
    return df

def plot_normalised_stacked_bar(data, labels, ticklabs, title="Stacked Bar Plot", colors=None, legend_title=None):
    """
    Plots a stacked vertical bar plot for given lists of data where each bar's total height is normalized to 1.

    Parameters:
        data (list of lists): Each sublist contains counts for elements, and each list represents a different category.
        labels (list): Labels for each of the elements (same length as each sublist in data).
        title (str): Title of the plot.
        colors (list): Colors for each element.
        legend_title (str): Title for the legend.

    Returns:
        A matplotlib bar plot.
    """
    data = np.array(data)
    data_normalized = data / data.sum(axis=1, keepdims=True)

    #if not colors:
        # Generate as many colors as there are labels if not provided
    #    colors = plt.cm.viridis(np.linspace(0, 1, len(labels)))

    # Setting up the plot
    fig, ax = plt.subplots(figsize = (5, 5))
    bar_width = 0.75
    indices = np.arange(data.shape[0])
    #print(indices)
    
    # Initialize the bottom array for the stacked bars
    bottom = np.zeros(data.shape[0])

    # Plot bars
    for i in range(data.shape[1]):
        ax.bar(indices, data_normalized[:, i], bar_width, bottom=bottom, color=colors[i], label=labels[i], linewidth = 1, edgecolor = "k")
        bottom += data_normalized[:, i]

    #ax.set_title("Where methods predict nothing")
    ax.set_xlabel('Method')
    ax.set_ylabel('p')
    ax.set_xticks(indices)
    ax.set_xticklabels(ticklabs, rotation = 45)
    ax.set_ylim(0, 1)
    plt.axhline(y = 0.5, linestyle = "--", color = "k")
    ax.legend(title=legend_title if legend_title else "Protein Group", bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.show()

def frequency_count(lst):
    return [lst.count("ELONGATED"), lst.count("ELONGATED+TINY"), lst.count("GLOBULAR"), lst.count("GLOBULAR+TINY")]

def plot_violin_with_swarm(df, col, col_lab, palette, my_order, FSIZE=(6,6), DPI=100, out = None, ymax = None, ymin = None, ymin_inc = 0, ymax_inc = 0):
    #Define a function to detect outliers
    def is_outlier(s):
        lower_limit = s.mean() - (s.std() * 4)
        upper_limit = s.mean() + (s.std() * 4)
        return ~s.between(lower_limit, upper_limit)
    
    # Apply the function to the dataframe
    df['Outlier'] = df.groupby('origin')[col].transform(is_outlier)

    plt.figure(figsize=FSIZE, dpi=DPI)
    
    # Create the violin plot with outliers removed
    sns.violinplot(x='origin', y=col, palette=palette, data=df[df['Outlier'] == False],
                   linewidth=0.5, linecolor="k", cut=0, inner='box', order=my_order, bw_adjust=0.7)

    # Overlay the swarm plot to show the outliers
    sns.swarmplot(x='origin', y=col, palette=palette, data=df[df['Outlier'] == True], order=my_order, size=2, edgecolor = "k", linewidth = 0.25)

    # Labeling the plot
    plt.xlabel('Method')
    plt.ylabel(col_lab)
    if ymax != None:
        if ymin == None:
            ymin = -0.5
        plt.ylim(ymin-ymin_inc,ymax+ymax_inc)

    else:
        ymax = df.query('Outlier == False')[col].max()
        if ymin == None:
            ymin = -0.5
        #print(ymax)
        plt.ylim(ymin-ymin_inc,ymax+ymax_inc)

    if out != None:
        plt.savefig(out)

    plt.show()


### VARIABLES ###

ASA_Wilke = {
    'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0,
    'CYS': 167.0, 'GLN': 225.0, 'GLU': 223.0, 'GLY': 104.0,
    'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0,
    'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0,
    'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0
}

palette = {
    'LIGYSIS': 'forestgreen',
    'VN-EGNN': 'orchid',
    'VN-EGNN-NR': 'orchid',
    'IF-SitePred': 'dodgerblue',
    'IF-SitePred-NR': 'dodgerblue',
    'IF-SitePred-rescored-NR': 'dodgerblue',
    'GrASP': 'firebrick',
    'PUResNet': 'royalblue',
    'PUResNet-AA': 'royalblue',
    'PUResNet+PRANK': 'royalblue',
    'DeepPocket-Segmented': 'sienna',
    'DeepPocket-Segmented-NR': 'sienna',
    'DeepPocket-Rescored': 'tan',
    'P2Rank+Cons': 'magenta',
    'P2Rank': 'orange',
    'fpocket+PRANK': 'gray',
    'fpocket': 'gray',
    'PocketFinder': 'blue',
    'PocketFinder-AA': 'blue',
    'PocketFinder-SS': 'blue',
    'PocketFinder+PRANK': 'blue',
    'Ligsite': 'red',
    'Ligsite-AA':  'red',
    'Ligsite-SS': 'red',
    'Ligsite+PRANK': 'red',
    'Surfnet': 'green',
    'Surfnet-AA': 'green',
    'Surfnet-SS': 'green',
    'Surfnet+PRANK': 'green',
 }

markers_dict = {
    "LIGYSIS": ".",
    "VN-EGNN": '8',
    "VN-EGNN-NR": '8',
    "IF-SitePred": 'D',
    "IF-SitePred-NR": 'D',
    "IF-SitePred-rescored-NR": 'D',
    "GrASP": 's',
    "PUResNet": 'p',
    "PUResNet-AA": 'p', 
    "PUResNet+PRANK": 'p',
    "DeepPocket-Segmented": 'H',
    "DeepPocket-Segmented-NR": 'H',
    "DeepPocket-Rescored": 'h',
    "P2Rank+Cons": 'v',
    "P2Rank": '^',
    "fpocket+PRANK": 'o',
    "fpocket": 'o',
    "PocketFinder": '*',
    "PocketFinder-AA": "*",
    "PocketFinder-SS": '*',
    "PocketFinder+PRANK": '*',
    "Ligsite": 'P',
    "Ligsite-AA": 'P',
    "Ligsite-SS": 'P',
    "Ligsite+PRANK": 'P',
    "Surfnet": 'X',
    "Surfnet-AA": 'X',
    'Surfnet-SS': 'X',
    "Surfnet+PRANK": 'X',
}


avg_max_sasa = statistics.mean(list(ASA_Wilke.values()))

rot_matrices = read_from_pickle("./results/PDB_rot_matrices.pkl")
orig_centroids = read_from_pickle("./results/PDB_orig_centroids.pkl")

clean_pdb_python_bin = "/Users/2394007/miniconda3/envs/ob/bin/python"
clean_pdb_bin = "./../PROGRAMS/clean_pdb.py"

