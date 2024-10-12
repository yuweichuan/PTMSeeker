# SeaPIC
SeaPIC is a screening tool that searches for post-translational modifications (PTMs) in the cross-linking mass spectrometry (XL-MS) data.
## Usage
### Environment
1. The Windows system is required for SeaPIC.
2. Fetch all the files in the repository.
3. Install Python (3.6 or above) and add it to the system path. (https://www.python.org/downloads/)
4. Install numpy, scipy, lxml, and pyteomics packages.
```bash
pip install numpy scipy lxml pyteomics
```
4. Test if the packages are successfully installed.
```bash
python
import numpy, scipy, lxml, pyteomics
```
you should observe something like below
```bash
C:\Users\zhouchen>python
Python 3.10.10 (tags/v3.10.10:aad5f6a, Feb  7 2023, 17:20:36) [MSC v.1929 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license" for more information.
>>> import numpy, scipy, lxml, pyteomics
>>>
```
### Quick start
Navigate to the SeaPIC directory and place your dataset and database in the designated path. You can start the data analysis by running the main.py with the parameter settings in params.txt. Upon completion, the program will generate CSV files containing the results for the task.
```bash
python main.py params.txt
```
The details of the input file params.txt are listed below.
- std_aa_mass
  - Standard amino acid masses list. Usually do not modify them if unnecessary.
- Parse_rule
  - Protein digestion rule. Users can choose "trypsin", "arg-c", etc.
- Fasta_path
  - Protein sequence database path.
- Data_path
  - MGF data path. Use zhouchen_xlms_ptm.ThermoXcalibur.opt to process the RAW data through the Mascot distiller.
- Max_length
  - Maximum peptide length to be digested.
- Min_length
  - Minimum peptide length to be digested.
- Miss_cleavage
  - Number of maximum allowed missed cleaves in the peptide.
- Num_max_mod
  - Maximum allowed number of known variable modifications in the peptides.
- Num_max_unknown
  - Maximum allowed number of unknown modifications in the peptides (maximum PTMs to be identified in one peptide).
- Link_site
  - Cross-linking reaction site. leave it empty to activate the linear peptide PTM identification task.
- xl_mass
  - Cross-linker residual mass. e.g. 138.068 for BS3.
- Known_fix_mod
  - Input known fix modification on peptides, such as C+57.02=160.03.
- Known_var_mod
  - Input known variable modification on peptides, such as M+15.99=147.04.
- upper_ptm_mass
  - Upper bound of PTMs mass to be considered.
- lower_ptm_mass
  - Lower bound of PTMs mass to be considered.
- upper_ptm_mass
  - Upper bound of PTMs mass to be considered.
- workers
  - The number of threads to use in the multiprocess.
- resolution
  - MS2 tolerance in Dalton. 0.01 is suggested. Usually do not exceed 0.02.

## Authors
czhouau@connect.ust.hk Chen Zhou

eeyu@ust.hk Weichuan Yu

## License
[MIT LICENSE]
