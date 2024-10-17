######################################################################################
# HAIRpred is developed for predicting and desigining the Human Antibody Binding Residues.  #
# It is developed by Prof G. P. S. Raghava's group.       #
# Please cite: https://webs.iiitd.edu.in/raghava/HAIRpred/                            #
######################################################################################

############### Packages to be Installed ##################
#pip install Bio == 1.7.1
#apt-get install dssp
#pip install gemmi == 0.6.7
#pip install scikit-learn==1.5.2
#pip install pandas == 2.1.4
#pip install numpy == 1.26.4
#pip install openpyxl == 3.1.5
#pip install transformers == 4.44.2
#pip install torch == 2.4.1
#pip install accelerate == 0.34.2


############## Importing necessary libraries ##################
import time
# Record start time
start_time = time.time()

import sys
import pandas as pd
import re
import os
import torch
import numpy as np
import shutil
import requests
import gemmi
from Bio.PDB.DSSP import DSSP  # Tool for calculating RSA
from Bio.PDB import PDBParser  # Parser for PDB structure files
import warnings
from urllib3.exceptions import InsecureRequestWarning
import copy
import pickle
import argparse
from sklearn.ensemble import RandomForestClassifier
from multiprocessing import Pool
from requests.exceptions import RequestException
import logging
import re
warnings.filterwarnings('ignore')
from transformers import logging
logging.set_verbosity_error()
from transformers.utils.logging import disable_progress_bar


################# Argument Parsing #####################

parser = argparse.ArgumentParser(description='Please provide following arguments') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence(s) in FASTA format or single sequence per line in single letter code")
parser.add_argument("-o", "--output",type=str, help="Output: File for saving results - It should have .xlsx as extension. Default : outfile.xlsx")
parser.add_argument("-j", "--job",type=int, choices = [1,2], help="Job Type: 1: Predict, 2: Design. Default : 1")
parser.add_argument("-m", "--model",type=int, choices = [1,2], help="(Only for Predict Module) Model Type: 1: RSA based RF, 2: RSA + PSSM ensemble model (Best Model). Default : 2")
parser.add_argument("-t","--threshold", type=float, help="Threshold: Value between 0 to 1. Default : 0.5")
parser.add_argument("-p","--path", type=str, help="Path for temporary folder")
args = parser.parse_args()


################# Functions #####################

# Function to read sequences from a FASTA file
def readseq(file):
    try:
        with open(file) as f:
            records = f.read()
        # Splitting the file content by '>' to process FASTA format
        records = records.split('>')[1:]
        seqid = []  # List for sequence IDs
        seq = []  # List for sequences
        special_chars_replaced = False  # Flag to check if any special characters were replaced
        # Process each sequence in the FASTA file
        for fasta in records:
            array = fasta.split('\n')
            # Extract the sequence name (ID) and clean up the sequence
            name = array[0]  # Keep the full sequence ID, even if it contains spaces
            original_name = name
            # Replace special characters with underscores
            name = re.sub(r'[ \(\)\|]', '_', name)
            if name != original_name:
                special_chars_replaced = True
            sequence = re.sub('[^ACDEFGHIKLMNPQRSTVWY-]', '', ''.join(array[1:]).upper())
            seqid.append(name)
            seq.append(sequence)
        
        # If no sequence IDs are found, handle as plain sequences line by line
        if len(seqid) == 0:
            with open(file, "r") as f:
                data1 = f.readlines()
            for each in data1:
                seq.append(each.replace('\n', ''))
            for i in range(1, len(seq) + 1):
                seqid.append("Seq_" + str(i))
        
        # Inform the user if special characters were replaced
        if special_chars_replaced:
            print("Note: Special characters (spaces, parentheses, '|') were found in sequence IDs. They have been replaced with underscores.")

        # Return DataFrame with sequence IDs and sequences
        df = pd.DataFrame({'seqid': seqid, 'seq': seq})
        return df
    
    # Handle file not found error
    except FileNotFoundError:
        print(f"Error: The file '{file}' was not found.")
        return None
    
    # Handle format errors or invalid data
    except (IndexError, ValueError) as e:
        print(f"Error: The input file '{file}' is in an incorrect format or contains invalid data.")
        return None

# Function to generate the PSSM file for a single sequence
def generate_pssm_for_sequence(seq_id, sequence, temp_dir, ncbi_dir):
    # Create a temporary FASTA file for the sequence
    temp_fasta_file = f"{temp_dir}/{seq_id}.fasta"
    with open(temp_fasta_file, "w") as temp_fasta:
        temp_fasta.write(f">{seq_id}\n{sequence}\n")
    # PSI-BLAST command to generate PSSM file
    cmd = (
        f"{ncbi_dir}/ncbi-blast-2.16.0+/bin/psiblast "
        f"-query {temp_fasta_file} -db {ncbi_dir}/swissprot/swissprot "
        f"-evalue 0.1 -word_size 3 -max_target_seqs 6000 -num_threads 10 "
        f"-gapopen 11 -gapextend 1 -matrix BLOSUM62 "
        f"-num_iterations 3 "
        f"-out_ascii_pssm {temp_dir}/pssm/{seq_id}.pssm"
        f" > /dev/null 2>&1"
    )
    os.system(cmd)  # Execute the command to generate the PSSM file
    os.remove(temp_fasta_file)  # Remove the temporary FASTA file

# Function to read the generated PSSM file
def get_pssm(pssm_id, temp_dir):
    # Read the PSSM file
    with open(f'{temp_dir}/pssm/{pssm_id}.pssm') as f:
        txt = f.read().splitlines()
        
        # Extract the PSSM matrix from the file
    pssm = []
    for i in range(3, len(txt) - 6):  # Skip header/footer lines
        ps = txt[i][10:-92].split()  # Extract relevant part of each line
        ps_int = [int(x) for x in ps]  # Convert to integers
        pssm.append(ps_int)
    return pssm

# Combined function to generate PSSM and fetch it
def generate_and_get_pssm(row, temp_dir, ncbi_dir):
    seq_id = row['seqid']  # Sequence ID
    sequence = row['seq']  # Sequence
    
    # Generate the PSSM file for the sequence
    generate_pssm_for_sequence(seq_id, sequence, temp_dir, ncbi_dir)
    
    # Get the generated PSSM data
    pssm_data = get_pssm(seq_id, temp_dir)
    
    return pssm_data

# Function to fetch PDB file using the ESMFold API
def fetch_pdb_file(seqid, sequence, save_path):
    # Suppress insecure request warnings
    warnings.simplefilter('ignore', InsecureRequestWarning)
    global job
    try:
        # Send a request to the ESMFold API with the sequence
        url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
        response = requests.post(url, data=sequence, timeout = 30, verify=False)
        response.raise_for_status()  # Raise an error for bad HTTP status codes
        if response.text:
            pdb_text = response.text  # Get the response text (PDB file content)
            # Save the PDB file to the specified path
            with open(save_path, 'w') as file:
                file.write(pdb_text)
            
            # Read the PDB file and save it in minimal PDB format
            st = gemmi.read_structure(save_path)
            st.write_minimal_pdb(save_path)
        else:
            print(f"Error: No response text received for seqid {seqid}")
    
    
    except RequestException as e:
        print(f"Request error for seqid {seqid}: {e}")
    except Exception as e:
            print(f"Error for seqid {seqid}: {e}")

# Function to calculate RSA using DSSP
def calculate_rsa(model, pdb_path):
    # Create a DSSP object for the PDB file
    try:
        dssp = DSSP(model, pdb_path, dssp='mkdssp')
        
        # Initialize an empty DataFrame for RSA values
        chain_residue_rsa = pd.DataFrame(columns=['Chain', 'Residue', 'RSA'])
        # Iterate over each residue and calculate RSA
        for (chain_id, residue_id) in dssp.keys():
            residue = dssp[(chain_id, residue_id)][1]
            if residue == 'X':  # Skip invalid residues
                continue
            rsa = dssp[(chain_id, residue_id)][3]  # Get RSA value
            chain_residue_rsa.loc[len(chain_residue_rsa)] = [chain_id, residue, rsa]
        
        return chain_residue_rsa
    except :
        if job ==2 :    # For design module - we have to take into account if DSSP for mutant is not possible
            chain_residue_rsa = pd.DataFrame(columns=['Chain', 'Residue', 'RSA'])
            return chain_residue_rsa
        if job ==1 : 
            print("Error with DSSP")

# Function to generate and fetch RSA values
def generate_and_get_rsa(row):
    seqid = row['seqid']  # Sequence ID
    pdb_path_unbound = f"./{seqid}_antigen.pdb"  # Path for the PDB file
    sequence = row['seq']  # Sequence
    if len(sequence)<=400:
       try : fetch_pdb_file(seqid, sequence, pdb_path_unbound)
       except : fetch_pdb_file_longer(seqid, sequence, pdb_path_unbound)
    else : fetch_pdb_file_longer(seqid, sequence, pdb_path_unbound)
    # Parse the PDB file
    p = PDBParser(QUIET=True)
    structure_unbound = p.get_structure("Antigen", pdb_path_unbound)
    model_unbound = structure_unbound[0]
    # Calculate RSA for the structure
    rsa_data_unbound = calculate_rsa(model_unbound, pdb_path_unbound)
    if len(rsa_data_unbound) == 0 and job ==2 : # For design module - we have to take into account if the mutant cannot be inputted to dssp
        return None
    rsa = rsa_data_unbound['RSA'].to_list()  # Extract RSA values
    os.remove(pdb_path_unbound)  # Remove the PDB file after use
    return rsa

import pandas as pd
import copy

def get_windows(row, model):     

    w = 7  # w = (window_length - 1) / 2; Best Model has windows = 15

    if model == 1:
        columns = ['Residue Number', 'Residue', 'RSA']
    elif model == 2:
        columns = ['Residue Number', 'Residue', 'RSA', 'PSSM']
    elif model == 3:
        columns = ['Residue Number', 'Residue', 'PSSM']
    res = pd.DataFrame(columns=columns)

    try:
        seq = row['seq']
        if model != 3:  # RSA is not used in model 3
            RSA = row['rsa']
            if len(seq) != len(RSA):
                raise ValueError("Length of sequence and RSA do not match")
        if model in [2, 3]:
            pssm = row['pssm']
        for i in range(len(seq)):
            residue_num = i + 1
            residue = seq[i]

            if model != 3:  # Initialize RSA window for models 1 and 2
                r = [0] * (2 * w + 1)

            if model in [2, 3]:  # Initialize PSSM window for models 2 and 3
                pssm_final = copy.deepcopy([[0] * 20] * (2 * w + 1))

            # Handle the case where i is less than w
            start_idx = max(0, i - w)
            end_idx = min(len(seq), i + w + 1)

            if model != 3:  # Get RSA slice for models 1 and 2
                r_slice = RSA[start_idx:end_idx]

            if model in [2, 3]:  # Get PSSM slice for models 2 and 3
                pssm_slice = pssm[start_idx:end_idx]

            # Determine the insertion point in the window
            insert_start = w - (i - start_idx)
            insert_end = insert_start + (end_idx - start_idx)

            # Insert the slices into the initialized windows
            if model != 3:
                r[insert_start:insert_end] = r_slice

            if model in [2, 3]:
                pssm_final[insert_start:insert_end] = pssm_slice

            # Store the result in the dataframe
            if model == 2:
                res.loc[len(res)] = [residue_num, residue, r, pssm_final]
            elif model == 3:
                res.loc[len(res)] = [residue_num, residue, pssm_final]
            else:
                res.loc[len(res)] = [residue_num, residue, r]
            

    except (KeyError, IndexError, ValueError) as e:
        print(f"Error processing: {e} for seqid {row['seqid']}")

    return res


#Function to run the RF models
def model_run(window_df, model, thres):
    global mod_rsa, mod_pssm
    # Initialize an empty dataframe to store results
    result_df = pd.DataFrame(columns=['Residue Number', 'Residue', 'HAIRpred Score', 'Prediction'])
    if model == 2: 
        for index, row in window_df.iterrows():
            # Get the RSA and PSSM for the residue
            rsa_window = row['RSA']  # Use the RSA window
            pssm_window = row['PSSM']  # Use the PSSM window (2D list)
            
            # Flatten the PSSM window for model input
            pssm_flat = np.array(pssm_window).flatten()

            # Ensure both models get valid input shapes
            rsa_input = np.array([rsa_window])
            pssm_input = np.array([pssm_flat])

            # Predict with both models
            rsa_prob = mod_rsa.predict_proba(rsa_input)[0,1]  # Assuming binary classification, we take the probability of class 1
            pssm_prob = mod_pssm.predict_proba(pssm_input)[0,1]  # Similarly for the PSSM model
            # Average the probabilities
            avg_prob = (rsa_prob + pssm_prob) / 2

            # Assign a label based on the threshold
            prediction = 1 if avg_prob >= thres else 0

            # Append the result to the dataframe
            result_df.loc[len(result_df)] = [row['Residue Number'], row['Residue'], avg_prob, prediction]

    if model == 1: 
         for index, row in window_df.iterrows():
            # Get the RSA and PSSM for the residue
            rsa_window = row['RSA']  # Use the RSA window

            # Ensure both models get valid input shapes
            rsa_input = np.array([rsa_window])

            # Predict with both models
            rsa_prob = mod_rsa.predict_proba(rsa_input)[0,1]  # Assuming binary classification, we take the probability of class 1

            # Assign a label based on the threshold
            prediction = 1 if rsa_prob >= thres else 0

            # Append the result to the dataframe
            result_df.loc[len(result_df)] = [row['Residue Number'], row['Residue'], rsa_prob, prediction]

    if model == 3: 
         for index, row in window_df.iterrows():
            # Get the RSA and PSSM for the residue
            pssm_window = row['PSSM']  # Use the RSA window

            # Flatten the PSSM window for model input
            pssm_flat = np.array(pssm_window).flatten()

            pssm_input = np.array([pssm_flat])

            # Predict with both models
            pssm_prob = mod_pssm.predict_proba(pssm_input)[0,1]  # Assuming binary classification, we take the probability of class 1

            # Assign a label based on the threshold
            prediction = 1 if pssm_prob >= thres else 0

            # Append the result to the dataframe
            result_df.loc[len(result_df)] = [row['Residue Number'], row['Residue'], pssm_prob, prediction]

    return result_df

# Predict Module to call different functions and give the result in an excel sheet
def predict(df, model, threshold, output):

    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        for i in range(len(df)):
            window_df = get_windows(df.loc[i], model)
            result_df = model_run(window_df, model, threshold)
            result_df.to_excel(writer, sheet_name=df.loc[i,'seqid'], index=False)

# Function for generating all possible mutants
def get_mutants(df):
    std = list("ACDEFGHIKLMNPQRSTVWY")  # Standard amino acids
    data = {'Seq_ID': [], 'seqid': [], 'seq': []}  # Initialize dictionary for storing data
    
    # Iterate through each row in the dataframe
    for k in range(len(df)):
        original_seq = df['seq'][k]  # Original sequence
        original_seqid = df['seqid'][k]  # Original sequence ID

        # Add original sequence to data
        data['Seq_ID'].append(original_seqid)
        data['seqid'].append(f'{original_seqid}_Original_Seq')
        data['seq'].append(original_seq)
        c=0
        # Generate mutants by replacing each residue with all other residues
        for i in range(len(original_seq)):
            for j in std:
                if original_seq[i] != j:  # Create mutant only if the residue differs
                    mutant_seq = original_seq[:i] + j + original_seq[i + 1:]  # Replace residue at position i with j
                    data['Seq_ID'].append(original_seqid)
                    data['seqid'].append(f'{original_seqid}_Mutant_' + str(c+1))
                    data['seq'].append(mutant_seq)
                    c=c+1
    
    # Create DataFrame directly from the collected data
    design_df = pd.DataFrame(data)
    
    # Return the final DataFrame containing all original and mutant sequences
    return design_df



####### ESMFold Functions (Accessed only when sequences are longer than 400 amino acids) ###########

# Function to load the model once
def load_esmfold_model():
    # Import necessary libraries
    from transformers import AutoTokenizer, EsmForProteinFolding
    disable_progress_bar()
    global esmfold_model, esmfold_tokenizer
    if esmfold_model is None:
        esmfold_tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
        esmfold_model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
        # Set model to use half precision to save memory
        esmfold_model.esm = esmfold_model.esm.half()
        torch.backends.cuda.matmul.allow_tf32 = True
        esmfold_model.trunk.set_chunk_size(64)
        if torch.cuda.is_available():
            esmfold_model = esmfold_model.cuda()


def fetch_pdb_file_longer(seqid, sequence, save_path):
    try:
        # Import necessary libraries
        from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
        from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

        global esmfold_model, esmfold_tokenizer
        if esmfold_model is None:
            load_esmfold_model()
        tokenized_input = esmfold_tokenizer(sequence, return_tensors="pt", add_special_tokens=False)['input_ids']

        # Move to GPU if available
        if torch.cuda.is_available():
            tokenized_input = tokenized_input.cuda()

        esmfold_model.eval()
        # Perform inference to get the model's outputs
        with torch.no_grad():
            outputs = esmfold_model(tokenized_input)

        # Convert atom14 positions to atom37 positions
        final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
        outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
        final_atom_positions = final_atom_positions.cpu().numpy()
        final_atom_mask = outputs["atom37_atom_exists"]

        # Generate PDB strings for each model output
        pdbs = []
        for i in range(outputs["aatype"].shape[0]):
            aa = outputs["aatype"][i]
            pred_pos = final_atom_positions[i]
            mask = final_atom_mask[i]
            resid = outputs["residue_index"][i] + 1
            pred = OFProtein(
                aatype=aa,
                atom_positions=pred_pos,
                atom_mask=mask,
                residue_index=resid,
                b_factors=outputs["plddt"][i],
                chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
            )
            pdbs.append(to_pdb(pred))

        # Write the PDB strings to the specified file
        with open(save_path, "w") as f:
            f.write("".join(pdbs))

        # Read the PDB file and save it in minimal PDB format
        st = gemmi.read_structure(save_path)
        st.write_minimal_pdb(save_path)

    except RequestException as e:
        print(f"Request error for seqid {seqid}: {e}")
    except Exception as e:
            print(f"Error in getting structure for seqid {seqid}: {e}")



################## Directory Paths ##########################
if args.path == None:
    temp_dir = "/Applications/XAMPP/xamppfiles/htdocs/hairpred/temp"
else: temp_dir = args.path # Directory to store temporary files
ncbi_dir = "/Applications/XAMPP/xamppfiles/htdocs/hairpred/pssm"  # Directory with NCBI BLAST and SwissProt database

temp_pssm_dir = f"{temp_dir}/pssm"
os.makedirs(temp_pssm_dir, exist_ok=True) # Ensure the PSSM directory exists

# Paths to the models in the GitHub repository
model_rsa_dir = "/Applications/XAMPP/xamppfiles/htdocs/hairpred/models/rf_bestmodel_rsa.pkl"
model_pssm_dir = "/Applications/XAMPP/xamppfiles/htdocs/hairpred/models/rf_bestmodel_pssm.pkl"

########## Initalizing ESMFold Model (only for sequences longer than 400 residues) #############

# Global variable to store the model, initialized to None
esmfold_model = None

################## Loading Models ###########################


# Load models once
with open(model_rsa_dir, 'rb') as f_rsa, open(model_pssm_dir, 'rb') as f_pssm:
    global mod_rsa, mod_pssm
    mod_rsa = pickle.load(f_rsa)
    mod_pssm = pickle.load(f_pssm)



################## Parameter initialization for command level arguments ########################

# Validate job argument for predict module
if args.job is not None and args.job not in [1, 2]:
    raise ValueError("Invalid value for job. In the predict module, the job must be 1 or 2.")

if args.output == None:
    output= "outfile.xlsx" 
else:
    output = args.output
         
# Threshold 
if args.threshold == None:
        threshold = 0.5
else:
        threshold= float(args.threshold)

# Job Type 
if args.job == None:
        job = int(1)
else:
        job = int(args.job)

# Model Type 
if args.model == None:
        model = int(2)
else:
        model = int(args.model)

# Validate that design module does not use job argument
if 'design' in sys.argv and args.job is not None:
    raise ValueError("The design module should not specify a job argument.")


if __name__ == '__main__':
    print('\n###############################################################################################')
    print('# Welcome to HAIRpred! #')
    print('# It is a Human Antibody Interation Prediction tool developed by Prof G. P. S. Raghava group. #')
    print('# Please cite: HAIRpred; available at https://webs.iiitd.edu.in/raghava/HAIRpred/  #')
    print('###############################################################################################\n')

    df = readseq(args.input)  # Read sequences from the FASTA file

    

    ##################### Prediction Module ########################

    if job == 1:
        print(f'\n======= Thanks for using Predict module of HAIRpred. Your results will be stored in file : {output} ===== \n')
        for i in range(len(df)):
            if len(df.loc[i,"seq"])>400 : 
                print("\nAtleast one of the sequences is longer than 400 residues. \nWe will be loading and running ESMFold on your device. It may take some time on devices without GPUs. \n")
                break
        
        if model == 1 :  #RSA only model
            try:
                """ Generate RSA for each sequence in parallel"""
                with Pool(processes=os.cpu_count()-1) as pool:
                    df['rsa'] = pool.starmap(generate_and_get_rsa, [(row,) for _, row in df.iterrows()])

                """  Prediction  """
                predict(df, model, threshold, output)

                print("\n=========Process Completed. Have a great day ahead! =============\n")
            except : print("RSA could not be generated for atleast one of the proteins. Please input foldable amino acid sequences. \n\n============ Have a great day ahead! ============= ")

        if model == 2 :  #RSA + PSSM ensemble model

            """ Generate RSA for each sequence in parallel"""
            try: 
                with Pool(processes=os.cpu_count()-1) as pool:
                    df['rsa'] = pool.starmap(generate_and_get_rsa, [(row,) for _, row in df.iterrows()])
            except : print("RSA could not be generated for atleast one of the proteins. Please input foldable amino acid sequences. \n\n============ Have a great day ahead! ============= ")

            try:
                """ Generate PSSM for each sequence"""
                # Run in parallel and assign PSSM data back to the dataframe
                with Pool(processes=os.cpu_count()-1) as pool:
                    df['pssm'] = pool.starmap(generate_and_get_pssm, [(row, temp_dir, ncbi_dir) for _, row in df.iterrows()])
                shutil.rmtree(temp_dir) # Remove the PSSM directory after use
                """  Prediction  """
                predict(df, model, threshold, output)
                print("\n=========Process Completed. Have a great day ahead! =============\n")
            except : print("PSSM could not be generated for atleast one of the proteins. Please select RSA based RF model. \n\n============ Have a great day ahead! =============")




    ##################### Design Module ########################

    if job == 2:
        try:
            print(f'\n======= Thanks for using Design module of HAIRpred. Your results will be stored in file : {output} ===== \n')
            df = get_mutants(df)
            """ Generate PSSM for each sequence"""
            # Run in parallel and assign PSSM data back to the dataframe
            with Pool(processes=os.cpu_count()-1) as pool:
                df['pssm'] = pool.starmap(generate_and_get_pssm, [(row, temp_dir, ncbi_dir) for _, row in df.iterrows()])
            shutil.rmtree(temp_dir) # Remove the PSSM directory after use
            """  Prediction  """
            predict(df, 3 , threshold, output)

            print("\n=========Process Completed. Have a great day ahead! =============\n")
        except: print("PSSM could not be generated for atleast one of the proteins. Please input a protein which has PSSM profile. \n\n============ Have a great day ahead! =============")





    # Record end time
    end_time = time.time()

    # Calculate elapsed time
    elapsed_time = end_time - start_time
    print(f"Time taken: {elapsed_time:.2f} seconds")

