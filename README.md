# HAIRpred : Prediction of human antibody interacting residues in an antigen from its primary structure



This repository contains the standalone code for HAIRpred prediction tool. HAIRpred is a human-specific antibody interacting residues predictor. To support the scientific community, we have developed a standalone software and web server for the tool. For a user-friendly access to HAIRpred, please use https://webs.iiitd.edu.in/raghava/hairpred/.

This github page was designed for users who have large proteins as inputs or large number of inputs. The source code available in this github page also gives users access to specific models, allowing for user customization.

## Reference
Sahni R, Kumar N, Raghava GPS. HAIRpred: Prediction of human antibody interacting residues in an antigen from its primary structure. Protein Sci. 2025 Aug;34(8):e70212.

## Installation (For Linux Users)

This Github repository has files which need Git-LFS for installation. To install Git-LFS : 

```bash
sudo apt-get install git-lfs
git lfs install
```

Then, clone this github repository : 

```bash
git clone https://github.com/Ruchir3003/HAIRpr.git
```

Next, we need to untar the tar.gz files :

```bash
tar -xvJf models.tar.xz
tar -xvfz pssm.tar.gz
```

### Install DSSP

To install DSSP, run the following command:

```bash 
apt-get install dssp
```
You can set up the environment using either `requirements.txt` (for pip users) or `environment.yml` (for Conda users).

### Using `requirements.txt`

```bash
pip install -r requirements.txt
```
### Using `environment.yml`

```bash
conda env create -f env.yml
```

## Installation (For Other Users)

You would need to install ncbi psi-blast files in the pssm and dssp for your system. 

### Install DSSP

You can install dssp by following https://github.com/cmbi/dssp

### Install PSSM

You can install system-specific ncbi psi-blast files from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

##  Usage
```
usage: python3 HAIRpred.py [-h]
                       [-i INPUT
                       [-o OUTPUT]
                       [-j {1,2}] 
                       [-m {1,2}] (Only for Predict Module)
                       [-t THRESHOLD]
```
```
Please provide following arguments for successful run

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or peptide sequence(s) in FASTA format
                        or single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -j {1,2}, --job {1,2}
                        Job Type: 1:Predict, 2: Design, by default 1
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.5
  -m {1,2}, --model {1,2}
                        Model Type: (Only for Predict Module) Model Type: 1: RSA based RF, 2: RSA + PSSM ensemble model (Best Model). Default : 2


**Input File:** It allow users to provide input in the FASTA format.

**Output File:** Program will save the results in the CSV format, in case user do not provide output file name, it will be stored in "outfile.csv".

**Threshold:** User should provide threshold between 0 and 1, by default its 0.5.

**Job:** User is allowed to choose between two different modules, such as, 1 for prediction and 2 for Designing , by default its 1.

**Model **: User is allowed to choose between two different models, such as, 1 for RSA based RF and 2 for RSA + PSSM ensemble RF, by default its 2.


