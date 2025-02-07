from flask import Flask, request, render_template_string
from dotenv import load_dotenv
import os
import re
from collections import Counter

# Load environment variables from .env file
load_dotenv()

# Create a Flask application instance
app = Flask(__name__)

# Configure the application with environment variables
app.config['SECRET_KEY'] = os.getenv('SECRET_KEY')
app.config['DATABASE_URL'] = os.getenv('DATABASE_URL')
app.config['API_KEY'] = os.getenv('API_KEY')

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
}

CODON_USAGE_TABLES = {
    'e_coli': {
        'ATA': 0.07, 'ATC': 0.42, 'ATT': 0.51, 'ATG': 1.00,  # Ile, Met
        'ACA': 0.13, 'ACC': 0.44, 'ACG': 0.27, 'ACT': 0.16,  # Thr
        'AAC': 0.55, 'AAT': 0.45, 'AAA': 0.74, 'AAG': 0.26,  # Asn, Lys
        'AGC': 0.28, 'AGT': 0.15, 'AGA': 0.07, 'AGG': 0.04,  # Ser, Arg
        'CTA': 0.04, 'CTC': 0.10, 'CTG': 0.50, 'CTT': 0.10,  # Leu
        'CCA': 0.19, 'CCC': 0.12, 'CCG': 0.52, 'CCT': 0.16,  # Pro
        'CAC': 0.57, 'CAT': 0.43, 'CAA': 0.34, 'CAG': 0.66,  # His, Gln
        'CGA': 0.06, 'CGC': 0.40, 'CGG': 0.10, 'CGT': 0.38,  # Arg
        'GTA': 0.15, 'GTC': 0.22, 'GTG': 0.37, 'GTT': 0.26,  # Val
        'GCA': 0.21, 'GCC': 0.27, 'GCG': 0.36, 'GCT': 0.16,  # Ala
        'GAC': 0.63, 'GAT': 0.37, 'GAA': 0.68, 'GAG': 0.32,  # Asp, Glu
        'GGA': 0.11, 'GGC': 0.41, 'GGG': 0.15, 'GGT': 0.33,  # Gly
        'TCA': 0.12, 'TCC': 0.15, 'TCG': 0.15, 'TCT': 0.15,  # Ser
        'TTC': 0.58, 'TTT': 0.42, 'TTA': 0.11, 'TTG': 0.13,  # Phe, Leu
        'TAC': 0.59, 'TAT': 0.41, 'TAA': 0.61, 'TAG': 0.09,  # Tyr, Stop
        'TGC': 0.56, 'TGT': 0.44, 'TGA': 0.30, 'TGG': 1.00   # Cys, Stop, Trp
    },
    'human': {
        'ATA': 0.16, 'ATC': 0.48, 'ATT': 0.36, 'ATG': 1.00,
        'ACA': 0.28, 'ACC': 0.36, 'ACG': 0.12, 'ACT': 0.24,
        'AAC': 0.54, 'AAT': 0.46, 'AAA': 0.42, 'AAG': 0.58,
        'AGC': 0.24, 'AGT': 0.15, 'AGA': 0.20, 'AGG': 0.20,
        'CTA': 0.07, 'CTC': 0.20, 'CTG': 0.41, 'CTT': 0.13,
        'CCA': 0.27, 'CCC': 0.33, 'CCG': 0.11, 'CCT': 0.29,
        'CAC': 0.58, 'CAT': 0.42, 'CAA': 0.25, 'CAG': 0.75,
        'CGA': 0.11, 'CGC': 0.19, 'CGG': 0.21, 'CGT': 0.08,
        'GTA': 0.11, 'GTC': 0.24, 'GTG': 0.47, 'GTT': 0.18,
        'GCA': 0.23, 'GCC': 0.40, 'GCG': 0.11, 'GCT': 0.26,
        'GAC': 0.54, 'GAT': 0.46, 'GAA': 0.42, 'GAG': 0.58,
        'GGA': 0.25, 'GGC': 0.34, 'GGG': 0.25, 'GGT': 0.16,
        'TCA': 0.15, 'TCC': 0.22, 'TCG': 0.05, 'TCT': 0.18,
        'TTC': 0.55, 'TTT': 0.45, 'TTA': 0.07, 'TTG': 0.13,
        'TAC': 0.56, 'TAT': 0.44, 'TAA': 0.28, 'TAG': 0.20,
        'TGC': 0.55, 'TGT': 0.45, 'TGA': 0.52, 'TGG': 1.00
    },
    'yeast': {
        'ATA': 0.27, 'ATC': 0.26, 'ATT': 0.47, 'ATG': 1.00,
        'ACA': 0.28, 'ACC': 0.21, 'ACG': 0.14, 'ACT': 0.37,
        'AAC': 0.43, 'AAT': 0.57, 'AAA': 0.58, 'AAG': 0.42,
        'AGC': 0.11, 'AGT': 0.16, 'AGA': 0.48, 'AGG': 0.21,
        'CTA': 0.14, 'CTC': 0.06, 'CTG': 0.11, 'CTT': 0.13,
        'CCA': 0.42, 'CCC': 0.15, 'CCG': 0.12, 'CCT': 0.31,
        'CAC': 0.36, 'CAT': 0.64, 'CAA': 0.69, 'CAG': 0.31,
        'CGA': 0.07, 'CGC': 0.06, 'CGG': 0.04, 'CGT': 0.15,
        'GTA': 0.21, 'GTC': 0.20, 'GTG': 0.19, 'GTT': 0.40,
        'GCA': 0.29, 'GCC': 0.21, 'GCG': 0.11, 'GCT': 0.39,
        'GAC': 0.37, 'GAT': 0.63, 'GAA': 0.70, 'GAG': 0.30,
        'GGA': 0.23, 'GGC': 0.19, 'GGG': 0.12, 'GGT': 0.46,
        'TCA': 0.21, 'TCC': 0.16, 'TCG': 0.10, 'TCT': 0.26,
        'TTC': 0.40, 'TTT': 0.60, 'TTA': 0.28, 'TTG': 0.28,
        'TAC': 0.43, 'TAT': 0.57, 'TAA': 0.47, 'TAG': 0.23,
        'TGC': 0.38, 'TGT': 0.62, 'TGA': 0.30, 'TGG': 1.00
    }
}

def parse_fasta(input_text):
    """Parse FASTA format text and return a dictionary of sequences."""
    sequences = {}
    current_header = ""
    current_sequence = []
    
    # Split the input into lines and process each line
    lines = [line.strip() for line in input_text.strip().split('\n')]
    
    for line in lines:
        if not line:  # Skip empty lines
            continue
        if line.startswith('>'):  # FASTA header line
            # If we were building a sequence, save it before starting a new one
            if current_header and current_sequence:
                sequences[current_header] = ''.join(current_sequence)
            # Start new sequence
            current_header = line[1:].strip()  # Remove '>' and whitespace
            current_sequence = []
        else:  # Sequence line
            # Remove any whitespace and convert to uppercase
            current_sequence.append(line.strip().upper())
    
    if current_header and current_sequence:
        sequences[current_header] = ''.join(current_sequence)
    
    return sequences

def process_input_sequence(input_text):
    input_text = input_text.strip()
    
    if input_text.startswith('>'):
        # This is FASTA format
        sequences = parse_fasta(input_text)
        if not sequences:
            raise ValueError("Invalid FASTA format")
        header, sequence = next(iter(sequences.items()))
        return sequence, header
    else:
        # This is a raw sequence
        sequence = ''.join(input_text.split())  # Remove all whitespace
        return sequence, None

# Kozak consensus sequence
kozak_regex = re.compile(r'(G|A)NN(A|G)TGATG')

def validate_dna(dna):
    """Validate the DNA sequence."""
    valid_nucleotides = set('ATCG')
    dna = ''.join(dna.split()).upper()
    return all(nucleotide in valid_nucleotides for nucleotide in dna)

def reverse_complement(dna):
    """Calculate the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(dna.upper()))

def calculate_sequence_complexity(dna):
    """Calculate sequence complexity using k-mer diversity."""
    k = 3  # using trinucleotides
    kmers = [dna[i:i+k] for i in range(len(dna)-k+1)]
    unique_kmers = len(set(kmers))
    max_possible = min(len(kmers), 4**k)  # 4^k is max possible k-mers for DNA
    return (unique_kmers / max_possible) * 100

def analyze_repeat_regions(dna):
    """Identify repeat regions in DNA sequence."""
    min_repeat_length = 4
    repeats = []
    
    for i in range(len(dna)-min_repeat_length):
        for j in range(min_repeat_length, min(20, len(dna)-i)):
            pattern = dna[i:i+j]
            if pattern in dna[i+j:]:
                repeats.append({
                    'pattern': pattern,
                    'length': len(pattern),
                    'position': i
                })
    
    # Remove overlapping repeats, keeping the longest ones
    filtered_repeats = []
    used_positions = set()
    
    for repeat in sorted(repeats, key=lambda x: -x['length']):
        pos = repeat['position']
        if not any(pos in range(p, p+r['length']) 
                  for p, r in zip(used_positions, filtered_repeats)):
            filtered_repeats.append(repeat)
            used_positions.add(pos)
            
    return filtered_repeats[:10]  # Return top 10 repeats

def find_orfs(dna):
    """Find all open reading frames in the given DNA sequence."""
    if len(dna) < 3:
        return []   
    pattern = re.compile(r'(?=(ATG(?:...)*?(?:TAA|TAG|TGA)))')
    return pattern.findall(dna)

def translate_dna(dna):
    """Translate a DNA sequence to a protein sequence."""
    protein = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        if codon in {'TAA', 'TAG', 'TGA'}:
            break
        amino_acid = gencode.get(codon, 'X')  # 'X' for unknown codons
        protein.append(amino_acid)
    return ''.join(protein)

def find_kozak_sequences(dna):
    """Find Kozak consensus sequences in the DNA."""
    return [match.start() for match in kozak_regex.finditer(dna)]

def calculate_cai(sequence: str) -> float:
        """Calculate Codon Adaptation Index."""
        codon_usage = {
            'ATA': 0.07, 'ATC': 0.25, 'ATT': 0.16, 'ATG': 1.0,
            'ACA': 0.14, 'ACC': 0.35, 'ACG': 0.09, 'ACT': 0.18,
            'AAC': 0.52, 'AAT': 0.48, 'AAA': 0.56, 'AAG': 0.44,
            'AGC': 0.31, 'AGT': 0.19, 'AGA': 0.15, 'AGG': 0.12,
            'CTA': 0.06, 'CTC': 0.22, 'CTG': 0.49, 'CTT': 0.12,
            'CCA': 0.20, 'CCC': 0.27, 'CCG': 0.10, 'CCT': 0.19,
            'CAC': 0.52, 'CAT': 0.48, 'CAA': 0.37, 'CAG': 0.63,
            'CGA': 0.06, 'CGC': 0.26, 'CGG': 0.12, 'CGT': 0.11,
            'GTA': 0.09, 'GTC': 0.22, 'GTG': 0.37, 'GTT': 0.17,
            'GCA': 0.17, 'GCC': 0.39, 'GCG': 0.08, 'GCT': 0.21,
            'GAC': 0.55, 'GAT': 0.45, 'GAA': 0.57, 'GAG': 0.43,
            'GGA': 0.19, 'GGC': 0.37, 'GGG': 0.13, 'GGT': 0.16,
            'TCA': 0.13, 'TCC': 0.22, 'TCG': 0.06, 'TCT': 0.19,
            'TTC': 0.52, 'TTT': 0.48, 'TTA': 0.08, 'TTG': 0.15,
            'TAC': 0.53, 'TAT': 0.47, 'TAA': 0.27, 'TAG': 0.04,
            'TGC': 0.52, 'TGT': 0.48, 'TGA': 0.68, 'TGG': 1.0
        }


        codons = [sequence[i:i+3] for i in range(0, len(sequence) - 2, 3)]
        cai_sum = 0
        for codon in codons:
            if codon in codon_usage:
                cai_sum += codon_usage[codon]
        return cai_sum / len(codons)

def predict_signal_peptide(protein):
    """
    Advanced signal peptide prediction using multiple features and position-specific scoring.
    
    Features considered:
    - N-region: Positive charges in first 5-8 residues
    - H-region: Hydrophobic core (8-12 residues)
    - C-region: Polar/small residues and cleavage site motifs
    """
    if len(protein) < 15:  # Minimum length for signal peptide
        return {"prediction": "No signal peptide detected", "confidence": 0, "details": "Sequence too short"}
        
    # Amino acid properties
    hydrophobic = set('AILMFWV')
    charged_positive = set('RK')
    charged_negative = set('DE')
    small_polar = set('GSTNQ')
    
    # Feature scores
    scores = {
        'n_region': 0,  # N-terminal charged region
        'h_region': 0,  # Hydrophobic core
        'c_region': 0,  # C-terminal region
        'cleavage': 0   # Cleavage site
    }
    
    # 1. N-region analysis (first 8 residues)
    n_region = protein[:8]
    n_positive_charges = sum(aa in charged_positive for aa in n_region)
    n_negative_charges = sum(aa in charged_negative for aa in n_region)
    scores['n_region'] = (n_positive_charges - n_negative_charges) * 0.5
    
    # 2. H-region analysis (residues 8-20)
    h_region = protein[8:20]
    hydrophobic_stretch = 0
    max_hydrophobic_stretch = 0
    
    for aa in h_region:
        if aa in hydrophobic:
            hydrophobic_stretch += 1
            max_hydrophobic_stretch = max(max_hydrophobic_stretch, hydrophobic_stretch)
        else:
            hydrophobic_stretch = 0
    
    scores['h_region'] = max_hydrophobic_stretch * 0.3
    
    # 3. C-region analysis (residues 20-30)
    if len(protein) >= 30:
        c_region = protein[20:30]
        small_polar_count = sum(aa in small_polar for aa in c_region)
        scores['c_region'] = small_polar_count * 0.2
        
        # 4. Cleavage site motif analysis
        # Common cleavage site motifs: A-X-A, G-X-A, etc.
        for i in range(len(c_region)-3):
            motif = c_region[i:i+3]
            if (motif[0] in 'AG' and motif[2] in 'A'):
                scores['cleavage'] += 1.0
    
    # Calculate total score and confidence
    total_score = sum(scores.values())
    max_possible_score = 5.0  # Theoretical maximum score
    confidence = (total_score / max_possible_score) * 100
    
    # Generate detailed analysis
    details = {
        'n_region_score': scores['n_region'],
        'h_region_score': scores['h_region'],
        'c_region_score': scores['c_region'],
        'cleavage_score': scores['cleavage'],
        'hydrophobic_core_length': max_hydrophobic_stretch,
        'positive_charges_n_term': n_positive_charges
    }
    
    # Make prediction
    if confidence >= 60:
        prediction = "Strong signal peptide detected"
    elif confidence >= 40:
        prediction = "Potential signal peptide detected"
    else:
        prediction = "No signal peptide detected"

    return {
        "prediction": prediction,
        "confidence": round(confidence, 2),
        "details": details
    }

def optimize_sequence(sequence, organism):
    """Optimize DNA sequence for expression in different organisms."""
    if organism not in CODON_USAGE_TABLES:
        raise ValueError(f"Optimization not available for {organism}")
    
    codon_table = CODON_USAGE_TABLES[organism]
    protein = translate_dna(sequence)
    optimized = ""
    
    # Get all possible codons for each amino acid
    aa_to_codons = {}
    for codon, aa in gencode.items():
        if aa not in aa_to_codons:
            aa_to_codons[aa] = []
        aa_to_codons[aa].append(codon)
    
    # For each amino acid in the protein sequence
    for aa in protein:
        if aa == '*':  # Stop codon
            possible_stops = ['TAA', 'TAG', 'TGA']
            best_stop = max(possible_stops, key=lambda x: codon_table.get(x, 0))
            optimized += best_stop
            break
            
        possible_codons = aa_to_codons.get(aa, [])
        if possible_codons:
            # Choose the codon with highest usage frequency
            best_codon = max(possible_codons, key=lambda x: codon_table.get(x, 0))
            optimized += best_codon
    
    # Calculate statistics
    original_cai = calculate_cai(sequence)
    optimized_cai = calculate_cai(optimized)
    gc_content_original = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    gc_content_optimized = (optimized.count('G') + optimized.count('C')) / len(optimized) * 100
    
    return {
        "optimized_sequence": optimized,
        "stats": {
            "original_cai": original_cai,
            "optimized_cai": optimized_cai,
            "improvement": f"{((optimized_cai - original_cai) / original_cai * 100):.1f}%",
            "gc_content_original": gc_content_original,
            "gc_content_optimized": gc_content_optimized
        }
    }

def analyze_dna(dna):
    """Analyze DNA sequence for ORFs, Kozak sequences, and translate to protein."""
    dna = ''.join(char for char in dna if not char.isspace())
    dna = dna.upper()
    
    if not validate_dna(dna):
        return {"error": "Invalid DNA sequence. Please use only A, T, C, and G."}
    
    """Check for continuous strings of one nucleotide"""
    if len(set(dna)) == 1:
        return {"error": "DNA sequence consists of a single nucleotide repeated."}

    complexity = calculate_sequence_complexity(dna)
    repeats = analyze_repeat_regions(dna)
    
    orfs = find_orfs(dna)
    if not orfs:
        return {"error": "No open reading frames found."}

    
    longest_orf = max(orfs, key=len)
    # optimization analysis for different organisms
    optimizations = {}
    for organism in ['e_coli', 'human', 'yeast']:
        try:
            optimizations[organism] = optimize_sequence(longest_orf, organism)
        except Exception as e:
            optimizations[organism] = {"error": str(e)}
    protein = translate_dna(longest_orf)
    kozak_positions = find_kozak_sequences(dna)
    cai = calculate_cai(longest_orf)
    signal_peptide_analysis = predict_signal_peptide(protein)
    gc_content = (dna.count('G') + dna.count('C')) / len(dna) * 100
    nucleotide_freq = {
        'A': dna.count('A'),
        'T': dna.count('T'),
        'G': dna.count('G'),
        'C': dna.count('C')
    }
    
    return {
        "longest_orf": longest_orf,
        "protein": protein,
        "kozak_positions": kozak_positions,
        "cai": cai,
        "signal_peptide": signal_peptide_analysis["prediction"],
        "signal_peptide_details": signal_peptide_analysis["details"],
        "signal_peptide_confidence": signal_peptide_analysis["confidence"],
        "gc_content": round(gc_content, 2),
        "nucleotide_frequency": nucleotide_freq,
        "sequence_length": len(dna),
        "reverse_complement": reverse_complement(dna),
        "optimizations": optimizations,
        "sequence_complexity": round(complexity, 2),
        "repeat_regions": repeats
    }

@app.route('/', methods=['GET', 'POST'])
def index():
    result = None
    error = None
    sequence_header = None
    if request.method == "POST":
        try:
            if 'dna_sequence' not in request.form:
                raise ValueError("No DNA sequence provided")
                
            # Get and sanitize input
            input_text = request.form['dna_sequence'].strip()
            
            # Validate input length
            if not input_text:
                raise ValueError("Sequence cannot be empty")
            
            # Process the sequence
            try:
                sequence, header = process_input_sequence(input_text)
                if header:
                    sequence_header = header
            except ValueError as e:
                raise ValueError(f"Invalid sequence format: {str(e)}") from e
            
            # basic XSS protection
            sequence = re.sub(r'[<>]', '', sequence)
            
            # Process the sequence
            result = analyze_dna(sequence)
            if sequence_header:
                result["sequence_header"] = sequence_header
            
            # Handle analysis errors
            if 'error' in result:
                error = result['error']
                result = None
                
        except ValueError as e:
            error = str(e)
        except Exception as e:
            error = "An unexpected error occurred. Please try again."
            app.logger.error(f"Error processing sequence: {str(e)}")


    return render_template_string('''   
        <!doctype html>
<html lang="en" class="dark">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>DNA2Protein Analysis Tool</title>
    <script src = "https://cdn.tailwindcss.com"></script>
    <link href="https://cdn.jsdelivr.net/npm/tailwindcss@2.2.19/dist/tailwind.min.css" rel="stylesheet">
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.2/css/all.min.css">
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <script>
        tailwind.config = {
            darkMode: 'class',
            theme: {
                extend: {
                    colors: {
                        primary: '#2C3E50',
                        secondary: '#34495E',
                        accent: '#576574'
                    },
                    animation: {
                        'dna-spin': 'spin 2s linear infinite',
                    }
                }
            }
        }
    </script>
    <style>
        @keyframes dna-spin {
            from { transform: rotate(0deg); }
            to { transform: rotate(360deg); }
        }
        
        .dna-spin {
            animation: dna-spin 2s linear infinite;
        }

        .result-value {
            scrollbar-width: thin;
            scrollbar-color: #4a5568 #1a202c;
        }
        
        .result-value::-webkit-scrollbar {
            width: 8px;
            height: 8px;
        }
        
        .result-value::-webkit-scrollbar-track {
            background: #1a202c;
        }
        
        .result-value::-webkit-scrollbar-thumb {
            background-color: #4a5568;
            border-radius: 4px;
        }

        #nucleotideChart {
            min-height: 300px;
            width: 100%;
    }

    </style>
</head>
<body class="bg-gradient-to-br from-gray-50 to-gray-100 dark:from-gray-900 dark:to-gray-800 min-h-screen text-gray-800 dark:text-gray-200">
    <div class="container mx-auto px-4 py-8 max-w-4xl">
        <h1 class="text-5xl font-extrabold text-center mb-12 bg-gradient-to-r from-primary to-secondary bg-clip-text text-transparent">DNA2PROTEIN</h1>
        
        <!-- Input Form -->
        <div class="bg-white dark:bg-gray-800 rounded-xl shadow-lg p-6 mb-8">
            <form method="post" id="dna-form" class="space-y-6">
                <div>
                    <label for="dna_sequence" class="block text-lg font-semibold mb-2">Enter DNA Sequence:</label>
                    <div class="mb-2 text-sm text-gray-600 dark:text-gray-400">
                        Accepts raw sequence or FASTA format (e.g., >sequence_name\nATGC...)
                    </div>
                    <textarea
                        id="dna_sequence" 
                        name="dna_sequence" 
                        required 
                        rows="6"
                        class="w-full p-4 bg-gray-50 dark:bg-gray-900 border border-gray-300 dark:border-gray-700 rounded-lg font-mono focus:ring-2 focus:ring-primary focus:border-transparent transition-all" 
                        placeholder=">sequence_name&#10;ATGCGATCGATCG"
                    ></textarea>
                </div>
                
                <div class="flex items-center gap-2 text-sm" id="sequence-validator">
                    <div class="w-3 h-3 rounded-full bg-gray-300"></div>
                    <span>Enter a valid DNA sequence</span>
                </div>

                <button type="submit" class="w-full bg-gradient-to-r from-primary to-secondary text-white font-semibold py-4 px-6 rounded-lg hover:opacity-90 transition-all flex items-center justify-center gap-3">
                    <i class="fas fa-dna hover:animate-dna-spin"></i>
                    Analyze
                </button>
            </form>
        </div>

        {% if error %}
        <div class="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded relative" role="alert">
            <strong class="font-bold">Error: </strong>
            <span class="block sm:inline">{{ error }}</span>
        </div>
        {% endif %}

        {% if result and not result.error %}
        <div class="bg-white dark:bg-gray-800 rounded-xl shadow-lg p-6 space-y-6">
            <h2 class="text-2xl font-bold border-b dark:border-gray-700 pb-4">Analysis Results</h2>
            
            {% if result.sequence_header %}
            <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                <div class="font-semibold mb-3 flex items-center gap-2">
                    <i class="fas fa-tag mr-2"></i>
                    <span>Sequence Name:</span>
                </div>
                <div class="font-mono bg-white dark:bg-gray-800 p-4 rounded-lg overflow-x-auto">
                    {{ result.sequence_header }}
                </div>
            </div>
            {% endif %}

            <!-- Result Items -->
            <div class="space-y-6">
                
                <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                    <div class="font-semibold mb-3 flex items-center gap-2">
                        <i class="fas fa-chart-pie mr-2"></i>
                        <span>Nucleotide Composition:</span>
                    </div>
                    <div style="height: 300px; position: relative;">
                        <canvas id="nucleotideChart"></canvas>
                    </div>
                </div>

                
            </div>

            <!-- Result Items -->
            <div class="space-y-6">
                <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                    <div class="font-semibold mb-3 flex items-center gap-2">
                        <i class="fas fa-dna mr-2 hover:animate-dna-spin"></i>
                        <span>Longest Open Reading Frame (ORF):</span>
                    </div>
                    <button onclick="copyToClipboard('{{ result.longest_orf }}', this)" 
                            class="text-gray-500 hover:text-gray-700 transition-colors">
                        <i class="fas fa-copy"></i>
                    </button>
                    <!-- Result Value -->
                    <div class="font-mono bg-white dark:bg-gray-800 p-4 rounded-lg overflow-x-auto">
                        {{ result.longest_orf }}
                    </div>
                </div>
            
            <!-- Result Items -->
            <div class="space-y-6">
                <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                    <div class="font-semibold mb-3 flex items-center gap-2">
                    <i class="fas fa-project-diagram mr-2"></i>
                    <span>Translated Protein:</span>
                </div>
                <button onclick="copyToClipboard('{{ result.protein }}', this)" 
                        class="text-gray-500 hover:text-gray-700 transition-colors">
                    <i class="fas fa-copy"></i>
                </button>
                <!-- Result Value -->
                <div class="font-mono bg-white dark:bg-gray-800 p-4 rounded-lg overflow-x-auto">{{ result.protein }}</div>
            </div>
            
            <!-- Result Items -->
            <div class="space-y-6">
                <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                    <div class="font-semibold mb-3 flex items-center gap-2">
                    <i class="fas fa-map-marker-alt mr-2"></i>
                    <span>Kozak Sequence Positions:</span>
                </div>

                <button onclick="copyToClipboard('{{ result.kozak_positions|join(', ') }}', this)" 
                        class="text-gray-500 hover:text-gray-700 transition-colors">
                    <i class="fas fa-copy"></i>
                </button>

                <!-- Result Value -->
                <div class="font-mono bg-white dark:bg-gray-800 p-4 rounded-lg overflow-x-auto">
                    {% if result.kozak_positions %}
                        {{ result.kozak_positions|join(', ') }}
                    {% else %}
                        No Kozak sequences found
                    {% endif %}
                </div>
            </div>
            
            <!-- Result Items -->
            <div class="space-y-6">
                <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                    <div class="font-semibold mb-3 flex items-center gap-2">
                    <i class="fas fa-chart-line mr-2"></i>
                    <span>Codon Adaptation Index (CAI):</span>
                </div>

                <button onclick="copyToClipboard('{{ result.cai }}', this)" 
                        class="text-gray-500 hover:text-gray-700 transition-colors">
                    <i class="fas fa-copy"></i>
                </button>

                <!-- Result Value -->
                <div class="font-mono bg-white dark:bg-gray-800 p-4 rounded-lg overflow-x-auto">{{ "%.2f"|format(result.cai) }}</div>
            </div>

            <!-- Result Items -->
            <div class="space-y-6">
                <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                    <div class="font-semibold mb-3 flex items-center gap-2">
                    <i class="fas fa-microscope mr-2"></i>
                    <span>Signal Peptide Prediction:</span>
                </div>

                <button onclick="copyToClipboard('{{ result.signal_peptide }}', this)" 
                        class="text-gray-500 hover:text-gray-700 transition-colors">
                    <i class="fas fa-copy"></i>
                </button>

                <!-- Result Value -->
                <div class="font-mono bg-white dark:bg-gray-800 p-4 rounded-lg overflow-x-auto">{{ result.signal_peptide }}</div>
            </div>

            <!-- Result Items -->
            <div class="space-y-6">
                <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                    <div class="font-semibold mb-3 flex items-center gap-2">
                    <i class="fas fa-percentage mr-2"></i>
                    <span>GC Content:</span>
                </div>

                <button onclick="copyToClipboard('{{ result.gc_content }}', this)" 
                        class="text-gray-500 hover:text-gray-700 transition-colors">
                    <i class="fas fa-copy"></i>
                </button>

                <!-- Result Value -->
                <div class="font-mono bg-white dark:bg-gray-800 p-4 rounded-lg overflow-x-auto">{{ result.gc_content }}%</div>
            </div>

            <!-- Result Items -->
            <div class="space-y-6">
                <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                    <div class="font-semibold mb-3 flex items-center gap-2">
                    <i class="fas fa-calculator mr-2"></i>
                    <span>Nucleotide Frequency:</span>
                </div>

                <button onclick="copyToClipboard('{{ result.nucleotide_frequency.A }},
                    {{ result.nucleotide_frequency.T }},
                    {{ result.nucleotide_frequency.G }},
                    {{ result.nucleotide_frequency.C }}', this)" 
                        class="text-gray-500 hover:text-gray-700 transition-colors">
                    <i class="fas fa-copy"></i>
                </button>

                <!-- Result Value -->
                <div class="font-mono bg-white dark:bg-gray-800 p-4 rounded-lg overflow-x-auto">
                    A: {{ result.nucleotide_frequency.A }},
                    T: {{ result.nucleotide_frequency.T }},
                    G: {{ result.nucleotide_frequency.G }},
                    C: {{ result.nucleotide_frequency.C }}
                </div>
            </div>

            <!-- Result Items -->
            <div class="space-y-6">
                <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                    <div class="font-semibold mb-3 flex items-center gap-2">
                    <i class="fas fa-ruler mr-2"></i>
                    </span>Sequence Length:</span>
                </div>

                <button onclick="copyToClipboard('{{ result.sequence_length }}', this)" 
                        class="text-gray-500 hover:text-gray-700 transition-colors">
                    <i class="fas fa-copy"></i>
                </button>

                <!-- Result Value -->
                <div class="font-mono bg-white dark:bg-gray-800 p-4 rounded-lg overflow-x-auto">{{ result.sequence_length }} bp</div>
            </div>

            <!-- Result Items -->
            <div class="space-y-6">
                <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                    <div class="font-semibold mb-3 flex items-center gap-2">
                    <i class="fas fa-arrow-left mr-2"></i>
                    <span>Reverse Complement:</span>
                </div>

                <button onclick="copyToClipboard('{{ result.reverse_complement }}', this)" 
                        class="text-gray-500 hover:text-gray-700 transition-colors">
                    <i class="fas fa-copy"></i>
                </button>

                <!-- Result Values -->
                <div class="font-mono bg-white dark:bg-gray-800 p-4 rounded-lg overflow-x-auto">{{ result.reverse_complement }}</div>
            </div>

            <!-- Result Items -->
            <div class="space-y-6">
                <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                    <div class="font-semibold mb-3 flex items-center gap-2">
                    <i class="fas fa-microscope mr-2"></i>
                    <span>Expression Optimization:</span>
                </div>

                

                <!-- Result Value -->
                <div class="font-mono bg-white dark:bg-gray-800 p-4 rounded-lg overflow-x-auto">
                    {% for organism, opt in result.optimizations.items() %}
                    <div class="mb-4">
                        <h4 class="font-semibold mb-2">{{ organism|title }} Optimization:</h4>
                        {% if opt.error %}
                            <p class="text-red-500">{{ opt.error }}</p>
                        {% else %}
                            <p>Original CAI: {{ "%.3f"|format(opt.stats.original_cai) }}</p>
                            <p>Optimized CAI: {{ "%.3f"|format(opt.stats.optimized_cai) }}</p>
                            <p>Improvement: {{ opt.stats.improvement }}</p>
                            <p>GC Content Change: {{ "%.1f"|format(opt.stats.gc_content_original) }}% → {{ "%.1f"|format(opt.stats.gc_content_optimized) }}%</p>
                            <details class="mt-2">
                                <summary class="cursor-pointer text-blue-500">View Optimized Sequence</summary>
                                <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                                <div class="font-semibold mb-3 flex items-center gap-2 overflow-x-auto">
                                    {{ opt.optimized_sequence }}
                                </div>
                                </div>
                            </details>
                        {% endif %}
                    </div>
                    {% endfor %}
                </div>
        </div>
        {% endif %}

        <footer class="mt-12 text-center text-gray-600 dark:text-gray-400">
            <p class="mb-2">Created by <a href="https://github.com/Bjorn99" class="hover:text-green-500 dark:text-white-400" target="_blank">Bjorn99</a></p>
            <p>
                <a href="https://github.com/Bjorn99/DNA2Protein" class="inline-flex items-center gap-2 text-white-800 hover:text-white-900 dark:text-white-400" target="_blank">
                    <i class="fab fa-github"></i>
                    View on GitHub
                </a>
            </p>
        </footer>
    </div>

    <script>
    document.addEventListener('DOMContentLoaded', function() {
        const input = document.getElementById('dna_sequence');
        const validator = document.getElementById('sequence-validator');
        const indicator = validator.querySelector('div');
        const text = validator.querySelector('span');

        function validateSequence(sequence) {
            if (!sequence) {
                return true;  // Empty is valid initially
            }

            // Handle FASTA format
            if (sequence.startsWith('>')) {
                // Split into lines and filter out empty lines
                const lines = sequence.split('\n').filter(line => line.trim());
                if (lines.length < 2) return true;  // Only header is present, still valid
                
                // Remove the header line and join the rest
                const dnaSequence = lines.slice(1).join('').replace(/\s+/g, '');
                return /^[ATCGatcg]+$/.test(dnaSequence);
            }

            // For regular sequence, just remove whitespace and check
            const cleanSequence = sequence.replace(/\s+/g, '');
            return /^[ATCGatcg]+$/.test(cleanSequence);
        }

        if (input && validator) {
            input.addEventListener('input', function() {
                const trimmedValue = this.value.trim();
                const isValid = validateSequence(trimmedValue);

                if (!trimmedValue) {
                    // Empty state
                    indicator.className = 'w-3 h-3 rounded-full bg-gray-300';
                    text.textContent = 'Enter a DNA sequence';
                    text.className = 'text-gray-600 dark:text-gray-400';
                } else {
                    // Valid or invalid state
                    if (isValid) {
                        indicator.className = 'w-3 h-3 rounded-full bg-green-500';
                        text.textContent = 'Valid sequence';
                        text.className = 'text-green-600 dark:text-green-400';
                    } else {
                        indicator.className = 'w-3 h-3 rounded-full bg-red-500';
                        text.textContent = 'Invalid sequence';
                        text.className = 'text-red-600 dark:text-red-400';
                    }
                }
            });
            
            // Initial validation state
            input.dispatchEvent(new Event('input'));
        }
    });
</script>

<script>
    // Separate chart initialization
    document.addEventListener('DOMContentLoaded', function() {
        {% if result and result.nucleotide_frequency %}
            const ctx = document.getElementById('nucleotideChart');
            if (ctx) {
                console.log('Creating chart with data:', {
                    A: {{ result.nucleotide_frequency.A }},
                    T: {{ result.nucleotide_frequency.T }},
                    G: {{ result.nucleotide_frequency.G }},
                    C: {{ result.nucleotide_frequency.C }}
                });
                
                new Chart(ctx, {
                    type: 'pie',
                    data: {
                        labels: ['A', 'T', 'G', 'C'],
                        datasets: [{
                            data: [
                                {{ result.nucleotide_frequency.A }},
                                {{ result.nucleotide_frequency.T }},
                                {{ result.nucleotide_frequency.G }},
                                {{ result.nucleotide_frequency.C }}
                            ],
                            backgroundColor: [
                                '#FF6384',
                                '#36A2EB',
                                '#FFCE56',
                                '#4BC0C0'
                            ]
                        }]
                    },
                    options: {
                        responsive: true,
                        maintainAspectRatio: false,
                        plugins: {
                            legend: {
                                position: 'top',
                            }
                        }
                    }
                });
            }
        {% endif %}
    });
</script>

<script>
function copyToClipboard(text, buttonElement) {
    navigator.clipboard.writeText(text).then(() => {
        const originalIcon = buttonElement.innerHTML;
        
        // Change to success icon
        buttonElement.innerHTML = '<i class="fas fa-check text-green-500"></i>';
        
        // Revert back after 2 seconds
        setTimeout(() => {
            buttonElement.innerHTML = originalIcon;
        }, 2000);
    }).catch(() => {
        // Error handling
        buttonElement.innerHTML = '<i class="fas fa-times text-red-500"></i>';
        setTimeout(() => {
            buttonElement.innerHTML = '<i class="fas fa-copy"></i>';
        }, 2000);
    });
}
</script>
</body>
</html>
    ''', result=result)         

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=False)