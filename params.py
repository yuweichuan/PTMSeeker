import json
std_aa_mass = {
    'G': 57.02146,
    'A': 71.03711,
    'S': 87.03203,
    'P': 97.05276,
    'V': 99.06841,
    'T': 101.04768,
    'C': 103.00919,
    # 'L': 113.08406,
    'I': 113.08406,
    'N': 114.04293,
    'D': 115.02694,
    'Q': 128.05858,
    'K': 128.09496,
    'E': 129.04259,
    'M': 131.04049,
    'H': 137.05891,
    'F': 147.06841,
    'R': 156.10111,
    'Y': 163.06333,
    'W': 186.07931
}

Parse_rule = 'trypsin'
Fasta_path = r"BSA.fasta"
Data_path = r"demo.mgf"
Max_length = 50
Min_length = 5
Miss_cleavage = 2
Num_max_mod = 2  # exclude fixed modifications
Num_max_unknown = 2  # max unknown ptm in one cross-linked peptide
Link_site = ['K']  # support multiple sites such as ['K', 'R']
xl_mass = 509.097  # CBDPS 509.097 138.068
Known_fix_mod = {'C': [103.00919 + 57.021464, 'c']}
Known_var_mod = {'M': [131.04049 + 15.9949, 'm']}
upper_ptm_mass = 250
lower_ptm_mass = -100
workers = 12  # multithreading workers
resolution = 0.01  # should better never be larger than 0.018 so the sake of distinguish ability


params = {'std_aa_mass': std_aa_mass, 'Parse_rule': Parse_rule, 'Fasta_path': Fasta_path, 'Data_path': Data_path,
          'Max_length': Max_length, 'Min_length': Min_length, 'Miss_cleavage': Miss_cleavage, 'Num_max_mod': Num_max_mod,
          'Num_max_unknown': Num_max_unknown, 'Link_site': Link_site, 'xl_mass': xl_mass, 'Known_fix_mod': Known_fix_mod,
          'Known_var_mod': Known_var_mod, 'upper_ptm_mass': upper_ptm_mass, 'lower_ptm_mass': lower_ptm_mass,
          'workers': workers, 'resolution': resolution}
json_object = json.dumps(params, indent=4)
with open('params.txt', 'w') as f:
    f.write(json_object)
