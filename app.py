from flask import Flask, request, render_template_string
from dotenv import load_dotenv
import os
import re
from collections import Counter
import requests
import json
from Bio import SeqIO, Entrez
from Bio.SeqUtils import GC, molecular_weight
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import joblib

# Load environment variables from .env file
load_dotenv()

# Create a Flask application instance
app = Flask(__name__)

# Configure the application with environment variables
app.config['SECRET_KEY'] = os.getenv('SECRET_KEY')
app.config['DATABASE_URL'] = os.getenv('DATABASE_URL')
app.config['API_KEY'] = os.getenv('83469fae42ffbc3846e3adbc8bc02fb09609')

# NCBI access
Entrez.api_key = os.getenv('83469fae42ffbc3846e3adbc8bc02fb09609')


# gencode = {
#     'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
#     'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
#     'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
#     'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
#     'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
#     'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
#     'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
#     'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
#     'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
#     'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
#     'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
#     'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
#     'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
#     'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
#     'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
#     'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
# }

# Complete genetic code dictionary with codon usage frequencies in E. coli
gencode = {
    # Isoleucine (I)
    'ATA': {'aa': 'I', 'freq': 0.07},
    'ATC': {'aa': 'I', 'freq': 0.48},
    'ATT': {'aa': 'I', 'freq': 0.45},
    
    # Methionine (M) - Start codon
    'ATG': {'aa': 'M', 'freq': 1.00},
    
    # Threonine (T)
    'ACA': {'aa': 'T', 'freq': 0.13},
    'ACC': {'aa': 'T', 'freq': 0.44},
    'ACG': {'aa': 'T', 'freq': 0.14},
    'ACT': {'aa': 'T', 'freq': 0.29},
    
    # Asparagine (N)
    'AAC': {'aa': 'N', 'freq': 0.51},
    'AAT': {'aa': 'N', 'freq': 0.49},
    
    # Lysine (K)
    'AAA': {'aa': 'K', 'freq': 0.74},
    'AAG': {'aa': 'K', 'freq': 0.26},
    
    # Serine (S)
    'AGC': {'aa': 'S', 'freq': 0.28},
    'AGT': {'aa': 'S', 'freq': 0.15},
    'TCA': {'aa': 'S', 'freq': 0.12},
    'TCC': {'aa': 'S', 'freq': 0.17},
    'TCG': {'aa': 'S', 'freq': 0.14},
    'TCT': {'aa': 'S', 'freq': 0.14},
    
    # Arginine (R)
    'AGA': {'aa': 'R', 'freq': 0.04},
    'AGG': {'aa': 'R', 'freq': 0.02},
    'CGA': {'aa': 'R', 'freq': 0.06},
    'CGC': {'aa': 'R', 'freq': 0.40},
    'CGG': {'aa': 'R', 'freq': 0.10},
    'CGT': {'aa': 'R', 'freq': 0.38},
    
    # Leucine (L)
    'CTA': {'aa': 'L', 'freq': 0.04},
    'CTC': {'aa': 'L', 'freq': 0.10},
    'CTG': {'aa': 'L', 'freq': 0.50},
    'CTT': {'aa': 'L', 'freq': 0.10},
    'TTA': {'aa': 'L', 'freq': 0.13},
    'TTG': {'aa': 'L', 'freq': 0.13},
    
    # Proline (P)
    'CCA': {'aa': 'P', 'freq': 0.19},
    'CCC': {'aa': 'P', 'freq': 0.12},
    'CCG': {'aa': 'P', 'freq': 0.52},
    'CCT': {'aa': 'P', 'freq': 0.17},
    
    # Histidine (H)
    'CAC': {'aa': 'H', 'freq': 0.43},
    'CAT': {'aa': 'H', 'freq': 0.57},
    
    # Glutamine (Q)
    'CAA': {'aa': 'Q', 'freq': 0.34},
    'CAG': {'aa': 'Q', 'freq': 0.66},
    
    # Valine (V)
    'GTA': {'aa': 'V', 'freq': 0.15},
    'GTC': {'aa': 'V', 'freq': 0.22},
    'GTG': {'aa': 'V', 'freq': 0.37},
    'GTT': {'aa': 'V', 'freq': 0.26},
    
    # Alanine (A)
    'GCA': {'aa': 'A', 'freq': 0.21},
    'GCC': {'aa': 'A', 'freq': 0.27},
    'GCG': {'aa': 'A', 'freq': 0.36},
    'GCT': {'aa': 'A', 'freq': 0.16},
    
    # Aspartic Acid (D)
    'GAC': {'aa': 'D', 'freq': 0.37},
    'GAT': {'aa': 'D', 'freq': 0.63},
    
    # Glutamic Acid (E)
    'GAA': {'aa': 'E', 'freq': 0.68},
    'GAG': {'aa': 'E', 'freq': 0.32},
    
    # Glycine (G)
    'GGA': {'aa': 'G', 'freq': 0.11},
    'GGC': {'aa': 'G', 'freq': 0.41},
    'GGG': {'aa': 'G', 'freq': 0.15},
    'GGT': {'aa': 'G', 'freq': 0.33},
    
    # Phenylalanine (F)
    'TTC': {'aa': 'F', 'freq': 0.42},
    'TTT': {'aa': 'F', 'freq': 0.58},
    
    # Tyrosine (Y)
    'TAC': {'aa': 'Y', 'freq': 0.41},
    'TAT': {'aa': 'Y', 'freq': 0.59},
    
    # Cysteine (C)
    'TGC': {'aa': 'C', 'freq': 0.56},
    'TGT': {'aa': 'C', 'freq': 0.44},
    
    # Tryptophan (W)
    'TGG': {'aa': 'W', 'freq': 1.00},
    
    # Stop codons
    'TAA': {'aa': '*', 'freq': 0.64},
    'TAG': {'aa': '*', 'freq': 0.07},
    'TGA': {'aa': '*', 'freq': 0.29}
}

class SignalPeptidePredictor:
    def __init__(self):
        self.model = self._initialize_model()
        self.amino_acid_properties = {
            'A': {'hydrophobicity': 1.8, 'charge': 0, 'size': 88.6},   # Alanine
            'R': {'hydrophobicity': -4.5, 'charge': 1, 'size': 173.4}, # Arginine
            'N': {'hydrophobicity': -3.5, 'charge': 0, 'size': 114.1}, # Asparagine
            'D': {'hydrophobicity': -3.5, 'charge': -1, 'size': 133.1}, # Aspartic Acid
            'C': {'hydrophobicity': 2.5, 'charge': 0, 'size': 121.2},   # Cysteine
            'E': {'hydrophobicity': -3.5, 'charge': -1, 'size': 147.1}, # Glutamic Acid
            'Q': {'hydrophobicity': -0.85, 'charge': 0, 'size': 146.2}, # Glutamine
            'G': {'hydrophobicity': -0.4, 'charge': 0, 'size': 75.1},   # Glycine
            'H': {'hydrophobicity': -3.2, 'charge': +1, 'size': 155.2}, # Histidine
            'I': {'hydrophobicity': 4.5, 'charge': 0, 'size': 131.2},   # Isoleucine
            'L': {'hydrophobicity': 3.8, 'charge': 0, 'size': 131.2},   # Leucine
            'K': {'hydrophobicity': -3.9, 'charge': +1, 'size': 146.2}, # Lysine
            'M': {'hydrophobicity': 1.9, 'charge': 0, 'size': 149.2},   # Methionine
            'F': {'hydrophobicity': 2.8, 'charge': 0, 'size': 165.2},   # Phenylalanine
            'P': {'hydrophobicity': -1.6, 'charge': 0, 'size': 115.1},   # Proline
            'S': {'hydrophobicity': -0.8, 'charge': 0, 'size': 105.1},   # Serine
            'T': {'hydrophobicity': -0.7, 'charge': 0, 'size': 119.1},   # Threonine
            'W': {'hydrophobicity': -0.9, 'charge': 0, 'size': 204.2},   # Tryptophan
            'Y': {'hydrophobicity': -1.3, 'charge': 0, 'size': 181.2},   # Tyrosine
            'V': {'hydrophobicity': 4.2, 'charge': 0, 'size': 117.1}    # Valine
        }

    def _initialize_model(self):
        try:
            return joblib.load('signal_peptide_model.pkl')
        except:
            # Return the basic model if the trained model is not available
            return RandomForestClassifier(n_estimators=100)

    def _extract_features(self, sequence):
        # Extract sequence features for ML model
        features = []
        window_size = 20

        # Sliding window analysis
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]
            hydrophobicity = sum(self.amino_acid_properties[aa]['hydrophobicity'] for aa in window)
            charge = sum(self.amino_acid_properties[aa]['charge'] for aa in window)
            size = sum(self.amino_acid_properties[aa]['size'] for aa in window)
            features.extend([hydrophobicity, charge, size])
            
        return np.array(features).reshape(1, -1)

    def predict(self, sequence):
        if len(sequence) < 30:
            return "Sequence too short for signal peptide prediction"
            
        features = self._extract_features(sequence[:30])
        prediction = self.model.predict_proba(features)[0]
        
        confidence = max(prediction) * 100
        is_signal = bool(prediction.argmax())

        return {
            'is_signal_peptide': is_signal,
            'confidence': confidence,
            'details': {
                'hydrophobicity_profile': self._calculate_hydrophobicity_profile(sequence[:30]),
                'charge_distribution': self._calculate_charge_distribution(sequence[:30])
            }
        }

    def _calculate_hydrophobicity_profile(self, sequence):
        return [self.amino_acid_properties[aa]['hydrophobicity'] for aa in sequence]

    def _calculate_charge_distribution(self, sequence):
        return [self.amino_acid_properties[aa]['charge'] for aa in sequence]

class BiologicalDatabaseIntegrator:
    def __init__(self):
        self.uniprot_api = "https://rest.uniprot.org/uniprotkb"
        self.pfam_api = "https://pfam.xfam.org/search/sequence"

    async def fetch_similar_sequences(self, sequence):
        try:
            # BLAST search through NCBI
            result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
            blast_records = BLAST.parse(result_handle)
            return [record for record in blast_records]
        except Exception as e:
            return f"Error in BLAST search: {str(e)}"

    async def get_protein_domains(self, sequence):
        try:
            response = requests.post(self.pfam_api, data={'seq': sequence})
            return response.json()
        except Exception as e:
            return f"Error in Pfam search: {str(e)}"

def analyze_dna(dna, include_extended_analysis=False):
    """Enhanced DNA sequence analysis function"""
    dna = ''.join(dna.split()).upper()
    
    if not validate_dna(dna):
        return {"error": "Invalid DNA sequence. Please use only A, T, C, and G."}
    
    try:
        # Basic analysis
        orfs = find_orfs(dna)
        if not orfs:
            return {"error": "No open reading frames found."}
        
        longest_orf = max(orfs, key=len)
        protein = translate_dna(longest_orf)
        
        # Enhanced analysis
        signal_peptide_predictor = SignalPeptidePredictor()
        signal_peptide_analysis = signal_peptide_predictor.predict(protein)
        
        # Calculate additional metrics
        gc_content = GC(dna)
        mol_weight = molecular_weight(dna)
        
        result = {
            "longest_orf": longest_orf,
            "protein": protein,
            "kozak_positions": find_kozak_sequences(dna),
            "cai": calculate_cai(longest_orf),
            "signal_peptide": signal_peptide_analysis,
            "gc_content": gc_content,
            "molecular_weight": mol_weight,
            "protein_length": len(protein)
        }

        if include_extended_analysis:
            db_integrator = BiologicalDatabaseIntegrator()
            result.update({
                "similar_sequences": await db_integrator.fetch_similar_sequences(protein),
                "protein_domains": await db_integrator.get_protein_domains(protein)
            })
        
        return result

    except Exception as e:
        return {"error": f"Analysis error: {str(e)}"}

# # Kozak consensus sequence
# kozak_regex = re.compile(r'(G|A)NN(A|G)TGATG')

# def validate_dna(dna):
#     """Validate the DNA sequence."""
#     valid_nucleotides = set('ATCG')
#     return all(nucleotide in valid_nucleotides for nucleotide in dna)

# def find_orfs(dna):
#     """Find all open reading frames in the given DNA sequence."""
#     if len(dna) < 3:
#         return []   
#     pattern = re.compile(r'(?=(ATG(?:...)*?(?:TAA|TAG|TGA)))')
#     return pattern.findall(dna)

# def translate_dna(dna):
#     """Translate a DNA sequence to a protein sequence."""
#     protein = []
#     for i in range(0, len(dna) - 2, 3):
#         codon = dna[i:i+3]
#         if codon in {'TAA', 'TAG', 'TGA'}:
#             break
#         amino_acid = gencode.get(codon, 'X')  # 'X' for unknown codons
#         protein.append(amino_acid)
#     return ''.join(protein)

# def find_kozak_sequences(dna):
#     """Find Kozak consensus sequences in the DNA."""
#     return [match.start() for match in kozak_regex.finditer(dna)]

# def calculate_cai(dna):
#     """Calculate the Codon Adaptation Index (CAI) for a DNA sequence."""
#     codon_counts = Counter(dna[i:i+3] for i in range(0, len(dna) - 2, 3))
#     total_codons = sum(codon_counts.values())
#     return sum(count / total_codons * len(gencode[codon]) for codon, count in codon_counts.items() if codon in gencode) / (len(dna) // 3)

# def predict_signal_peptide(protein):
#     """Simple prediction of signal peptide presence based on N-terminal sequence."""
#     n_terminal = protein[:30]  # Consider first 30 amino acids
#     if n_terminal.count('L') + n_terminal.count('A') + n_terminal.count('V') > 10:
#         return "Potential signal peptide detected"
#     return "No signal peptide detected"

# def analyze_dna(dna):
#     """Analyze DNA sequence for ORFs, Kozak sequences, and translate to protein."""
#     dna = ''.join(dna.split()).upper()
    
#     if not validate_dna(dna):
#         return {"error": "Invalid DNA sequence. Please use only A, T, C, and G."}
    
#     """Check for continuous strings of one nucleotide"""
#     if len(set(dna)) == 1:
#         return {"error": "DNA sequence consists of a single nucleotide repeated."}
    
#     orfs = find_orfs(dna)
#     if not orfs:
#         return {"error": "No open reading frames found."}
    
#     longest_orf = max(orfs, key=len)
#     protein = translate_dna(longest_orf)
#     kozak_positions = find_kozak_sequences(dna)
#     cai = calculate_cai(longest_orf)
#     signal_peptide = predict_signal_peptide(protein)
    
#     return {
#         "longest_orf": longest_orf,
#         "protein": protein,
#         "kozak_positions": kozak_positions,
#         "cai": cai,
#         "signal_peptide": signal_peptide
#     }

@app.route('/', methods=['GET', 'POST'])
def index():
    result = None
    if request.method == "POST":
        dna_sequence = request.form['dna_sequence']
        include_extended = request.form.get('include_extended', False)
        result = analyze_dna(dna_sequence, include_extended)

    return render_template_string('''   
        <!doctype html>
        <html lang="en">
        <head>
            <meta charset="utf-8">
            <meta name="viewport" content="width=device-width, initial-scale=1">
            <title>DNA2Protein</title>
            <link href="https://cdn.jsdelivr.net/npm/tailwindcss@2.2.19/dist/tailwind.min.css" rel="stylesheet">
            <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0-beta3/css/all.min.css">
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            <style>
                body {
                    background-color: #f0f4f8;
                    color: #2d3748;
                    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                }
                .container {
                    max-width: 800px;
                    margin: 0 auto;
                    padding: 2rem;
                }
                .title {
                    color: #2b6cb0;
                    font-size: 2.5rem;
                    font-weight: bold;
                    text-align: center;
                    margin-bottom: 2rem;
                    text-shadow: 1px 1px 2px rgba(0,0,0,0.1);
                }
                .input-form {
                    background-color: white;
                    padding: 2rem;
                    border-radius: 12px;
                    box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05);
                    display: flex;
                    flex-direction: column;
                    align-items: center;
                }
                .input-label {
                    font-weight: bold;
                    margin-bottom: 0.5rem;
                    color: #4a5568;
                    font-size: 1.1rem;
                }
                .input-field {
                    width: 100%;
                    padding: 0.75rem;
                    border: 2px solid #cbd5e0;
                    border-radius: 8px;
                    font-size: 1rem;
                    transition: all 0.3s ease;
                    margin-bottom: 1.5rem;
                }
                .input-field:focus {
                    border-color: #4299e1;
                    outline: none;
                    box-shadow: 0 0 0 3px rgba(66, 153, 225, 0.5);
                }
                .submit-button {
                    position: relative;
                    overflow: hidden;
                    transition: all 0.3s ease;
                    background-color: #4299e1;
                    color: white;
                    font-weight: bold;
                    padding: 0.75rem 2rem;
                    border-radius: 8px;
                    font-size: 1.1rem;
                    text-transform: uppercase;
                    letter-spacing: 0.05em;
                }
                .submit-button:hover {
                    background-color: #3182ce;
                    transform: translateY(-2px);
                    box-shadow: 0 4px 6px rgba(66, 153, 225, 0.5);
                }
                .submit-button::before {
                    content: '';
                    position: absolute;
                    top: -50%;
                    left: -50%;
                    width: 200%;
                    height: 200%;
                    background: radial-gradient(circle, rgba(255,255,255,0.3) 0%, rgba(255,255,255,0) 70%);
                    transform: scale(0);
                    transition: transform 0.6s ease-out;
                }
                .submit-button:hover::before {
                    transform: scale(1);
                }
                .submit-button i {
                    display: inline-block;
                    margin-right: 0.5rem;
                }
                .submit-button:hover i {
                    animation: dna-spin 2s linear infinite;
                }
                .submit-button:active {
                    transform: translateY(0);
                    animation: button-pulse 0.3s ease-out;
                }
                @keyframes dna-spin {
                    0% { transform: rotate(0deg); }
                    100% { transform: rotate(360deg); }
                }
                @keyframes button-pulse {
                    0%, 100% { transform: scale(1); }
                    50% { transform: scale(1.05); }
                }
                .result-box {
                    background-color: white;
                    border-radius: 12px;
                    padding: 1.5rem;
                    margin-top: 2rem;
                    box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05);
                    transition: all 0.3s ease;
                }
                .result-box:hover {
                    transform: translateY(-5px);
                    box-shadow: 0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04);
                }
                .result-title {
                    font-size: 1.5rem;
                    font-weight: bold;
                    color: #2b6cb0;
                    margin-bottom: 1rem;
                    border-bottom: 2px solid #e2e8f0;
                    padding-bottom: 0.5rem;
                }
                .result-item {
                    margin-bottom: 1rem;
                    transition: all 0.3s ease;
                }
                .result-item:hover {
                    transform: translateX(5px);
                }
                .result-label {
                    font-weight: bold;
                    color: #4a5568;
                    margin-bottom: 0.25rem;
                }
                .result-value {
                    background-color: #edf2f7;
                    padding: 0.75rem;
                    border-radius: 8px;
                    font-family: 'Courier New', Courier, monospace;
                    word-break: break-all;
                    font-size: 0.9rem;
                    transition: all 0.3s ease;
                }
                .result-value:hover {
                    background-color: #e2e8f0;
                }
                .error-message {
                    color: #e53e3e;
                    font-weight: bold;
                    margin-top: 1rem;
                    padding: 1rem;
                    border-radius: 8px;
                    background-color: #fff5f5;
                    border: 1px solid #feb2b2;
                }
            </style><style>
                body {
                    background-color: #f0f4f8;
                    color: #2d3748;
                    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                }
                .container {
                    max-width: 800px;
                    margin: 0 auto;
                    padding: 2rem;
                }
                .title {
                    color: #2b6cb0;
                    font-size: 2.5rem;
                    font-weight: bold;
                    text-align: center;
                    margin-bottom: 2rem;
                    text-shadow: 1px 1px 2px rgba(0,0,0,0.1);
                }
                .input-form {
                    background-color: white;
                    padding: 2rem;
                    border-radius: 12px;
                    box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05);
                    display: flex;
                    flex-direction: column;
                    align-items: center;
                }
                .input-label {
                    font-weight: bold;
                    margin-bottom: 0.5rem;
                    color: #4a5568;
                    font-size: 1.1rem;
                }
                .input-field {
                    width: 100%;
                    padding: 0.75rem;
                    border: 2px solid #cbd5e0;
                    border-radius: 8px;
                    font-size: 1rem;
                    transition: all 0.3s ease;
                    margin-bottom: 1.5rem;
                }
                .input-field:focus {
                    border-color: #4299e1;
                    outline: none;
                    box-shadow: 0 0 0 3px rgba(66, 153, 225, 0.5);
                }
                .submit-button {
                    position: relative;
                    overflow: hidden;
                    transition: all 0.3s ease;
                    background-color: #4299e1;
                    color: white;
                    font-weight: bold;
                    padding: 0.75rem 2rem;
                    border-radius: 8px;
                    font-size: 1.1rem;
                    text-transform: uppercase;
                    letter-spacing: 0.05em;
                }
                .submit-button:hover {
                    background-color: #3182ce;
                    transform: translateY(-2px);
                    box-shadow: 0 4px 6px rgba(66, 153, 225, 0.5);
                }
                .submit-button::before {
                    content: '';
                    position: absolute;
                    top: -50%;
                    left: -50%;
                    width: 200%;
                    height: 200%;
                    background: radial-gradient(circle, rgba(255,255,255,0.3) 0%, rgba(255,255,255,0) 70%);
                    transform: scale(0);
                    transition: transform 0.6s ease-out;
                }
                .submit-button:hover::before {
                    transform: scale(1);
                }
                .submit-button i {
                    display: inline-block;
                    margin-right: 0.5rem;
                }
                .submit-button:hover i {
                    animation: dna-spin 2s linear infinite;
                }
                .submit-button:active {
                    transform: translateY(0);
                    animation: button-pulse 0.3s ease-out;
                }
                @keyframes dna-spin {
                    0% { transform: rotate(0deg); }
                    100% { transform: rotate(360deg); }
                }
                @keyframes button-pulse {
                    0%, 100% { transform: scale(1); }
                    50% { transform: scale(1.05); }
                }
                .result-box {
                    background-color: white;
                    border-radius: 12px;
                    padding: 1.5rem;
                    margin-top: 2rem;
                    box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05);
                    transition: all 0.3s ease;
                }
                .result-box:hover {
                    transform: translateY(-5px);
                    box-shadow: 0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04);
                }
                .result-title {
                    font-size: 1.5rem;
                    font-weight: bold;
                    color: #2b6cb0;
                    margin-bottom: 1rem;
                    border-bottom: 2px solid #e2e8f0;
                    padding-bottom: 0.5rem;
                }
                .result-item {
                    margin-bottom: 1rem;
                    transition: all 0.3s ease;
                }
                .result-item:hover {
                    transform: translateX(5px);
                }
                .result-label {
                    font-weight: bold;
                    color: #4a5568;
                    margin-bottom: 0.25rem;
                }
                .result-value {
                    background-color: #edf2f7;
                    padding: 0.75rem;
                    border-radius: 8px;
                    font-family: 'Courier New', Courier, monospace;
                    word-break: break-all;
                    font-size: 0.9rem;
                    transition: all 0.3s ease;
                }
                .result-value:hover {
                    background-color: #e2e8f0;
                }
                .error-message {
                    color: #e53e3e;
                    font-weight: bold;
                    margin-top: 1rem;
                    padding: 1rem;
                    border-radius: 8px;
                    background-color: #fff5f5;
                    border: 1px solid #feb2b2;
                }
                .footer {
                    margin-top: 2rem;
                    text-align: center;
                    font-size: 0.9rem;
                    color: #4a5568;
                }
                .footer a {
                    color: #4299e1;
                    text-decoration: none;
                    transition: color 0.3s ease;
                }
                .footer a:hover {
                    color: #2b6cb0;
                    text-decoration: underline;
                }
                .feature-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
                gap: 1.5rem;
                margin-top: 2rem;
            }
            
            .feature-card {
                background: white;
                border-radius: 12px;
                padding: 1.5rem;
                box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
                transition: transform 0.3s ease;
            }
            
            .feature-card:hover {
                transform: translateY(-5px);
            }
            
            .visualization {
                width: 100%;
                height: 400px;
                margin-top: 2rem;
            }
            
            .loading-overlay {
                position: fixed;
                top: 0;
                left: 0;
                width: 100%;
                height: 100%;
                background: rgba(255, 255, 255, 0.9);
                display: flex;
                justify-content: center;
                align-items: center;
                z-index: 1000;
            }
            
            .loading-spinner {
                border: 4px solid #f3f3f3;
                border-top: 4px solid #3498db;
                border-radius: 50%;
                width: 50px;
                height: 50px;
                animation: spin 1s linear infinite;
            }
            
            @keyframes spin {
                0% { transform: rotate(0deg); }
                100% { transform: rotate(360deg); }
            }
            </style>
        </head>
        <body class="bg-gradient-to-br from-blue-50 to-indigo-100 min-h-screen">
            <div class="container">
                <h1 class="title">DNA2Protein</h1>
                <form method="post" class="input-form">
                    <label for="dna_sequence" class="input-label">Enter DNA Sequence:</label>
                    <input type="text" id="dna_sequence" name="dna_sequence" required class="input-field" placeholder="e.g., ATGCGATCGATCG"     >
                    <button type="submit" class="submit-button bg-blue-400 hover:bg-green-600 text-white font-bold py-2 px-4 rounded transition duration-300 ease-in-out transform hover:-translate-y-1 hover:shadow-lg">
                            <i class="fas fa-dna mr-2"></i> Analyze DNA
                    </button>
                </form>
                
                {% if result %}
                    <div class="result-box">
                        <h2 class="result-title">Analysis Result</h2>
                        {% if result.error %}
                            <div class="error-message">{{ result.error }}</div>
                        {% else %}
                            <div class="result-item">
                                <div class="result-label">Longest Open Reading Frame (ORF):</div>
                                <div class="result-value">{{ result.longest_orf }}</div>
                            </div>
                            <div class="result-item">
                                <div class="result-label">Translated Protein:</div>
                                <div class="result-value">{{ result.protein }}</div>
                            </div>
                            <div class="result-item">
                                <div class="result-label">Kozak Sequence Positions:</div>
                                <div class="result-value">
                                    {% if result.kozak_positions %}
                                        {{ result.kozak_positions|join(', ') }}
                                    {% else %}
                                        No Kozak sequences found
                                    {% endif %}
                                </div>
                            </div>
                            <div class="result-item">
                                <div class="result-label">Codon Adaptation Index (CAI):</div>
                                <div class="result-value">{{ "%.2f"|format(result.cai) }}</div>
                            </div>
                            <div class="result-item">
                                <div class="result-label">Signal Peptide Prediction:</div>
                                <div class="result-value">{{ result.signal_peptide }}</div>
                            </div>
                        {% endif %}

                        <div class="footer">
                            <p>Created by <a href="https://github.com/Bjorn99" target="_blank">Bjorn99</a></p>
                            <p>Have issues or want to contribute? Visit the <a href="https://github.com/Bjorn99/DNA2Protein" target="_blank">GitHub repository</a>.</p>
                        </div>
                    </div>
                {% endif %}
            </div>
        </body>
        </html>
    ''', result=result)         

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=True)