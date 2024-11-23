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

# Kozak consensus sequence
kozak_regex = re.compile(r'(G|A)NN(A|G)TGATG')

def validate_dna(dna):
    """Validate the DNA sequence."""
    valid_nucleotides = set('ATCG')
    return all(nucleotide in valid_nucleotides for nucleotide in dna)

def reverse_complement(dna):
    """Calculate the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(dna.upper()))

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
    """Simple prediction of signal peptide presence based on N-terminal sequence."""
    n_terminal = protein[:30]  # Consider first 30 amino acids
    if n_terminal.count('L') + n_terminal.count('A') + n_terminal.count('V') > 10:
        return "Potential signal peptide detected"
    return "No signal peptide detected"

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
    signal_peptide = predict_signal_peptide(protein)
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
        "signal_peptide": signal_peptide,
        "gc_content": round(gc_content, 2),
        "nucleotide_frequency": nucleotide_freq,
        "sequence_length": len(dna),
        "reverse_complement": reverse_complement(dna),
        "optimizations": optimizations
    }

@app.route('/', methods=['GET', 'POST'])
def index():
    result = None
    error = None
    if request.method == "POST":
        try:
            # Check if form data exists
            if 'dna_sequence' not in request.form:
                raise ValueError("No DNA sequence provided")
                
            # Get and sanitize input
            dna_sequence = request.form['dna_sequence'].strip()
            
            # Validate input length
            if not dna_sequence:
                raise ValueError("DNA sequence cannot be empty")
            
            # if len(dna_sequence) > 10000:  # Adjust limit as needed
            #     raise ValueError("DNA sequence is too long (maximum 10000 bases)")
                
            # Add basic XSS protection
            dna_sequence = re.sub(r'[<>]', '', dna_sequence)
            
            # Process the sequence
            result = analyze_dna(dna_sequence)
            
            # Handle analysis errors
            if 'error' in result:
                error = result['error']
                result = None
                
        except ValueError as e:
            error = str(e)
        except Exception as e:
            error = "An unexpected error occurred. Please try again."
            # Log the actual error for debugging
            app.logger.error(f"Error processing DNA sequence: {str(e)}")

    return render_template_string('''   
        <!doctype html>
<html lang="en" class="dark">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>DNA2Protein Analysis Tool</title>
    <script src = "https://cdn.tailwindcss.com"></script>
    <link href="https://cdn.jsdelivr.net/npm/tailwindcss@2.2.19/dist/tailwind.min.css" rel="stylesheet">
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.7.1/css/all.min.css">
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
                    <input type="text" 
                        id="dna_sequence" 
                        name="dna_sequence" 
                        required 
                        class="w-full p-4 bg-gray-50 dark:bg-gray-900 border border-gray-300 dark:border-gray-700 rounded-lg font-mono focus:ring-2 focus:ring-primary focus:border-transparent transition-all" 
                        placeholder="e.g., ATGCGATCGATCG"
                        title="Please enter valid DNA sequence (A, T, C, G only)">
                </div>
                
                
                <div class="flex items-center gap-2 text-sm" id="sequence-validator">
                    <div class="w-3 h-3 rounded-full bg-gray-300"></div>
                    <span class="validator-text">Enter a valid DNA sequence</span>
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
            
            {% if result.error %}
            <div class="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded relative" role="alert">
                <strong class="font-bold">Error: </strong>
                <span class="block sm:inline">{{ result.error }}</span>
            </div>
            {% else %}
            <!-- Result Items -->
            <div class="space-y-6">
                <div class="result-item bg-gray-50 dark:bg-gray-900 rounded-lg p-6">
                    <div class="font-semibold mb-3 flex items-center gap-2">
                        <i class="fas fa-dna mr-2 hover:animate-dna-spin"></i>
                        <span>Longest Open Reading Frame (ORF):</span>
                    </div>
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
            const text = document.querySelector('span');

            function validateSequence(sequence) {
                return /^[ATCGatcg\\s]+$/.test(sequence)
            }

            input.addEventListener('input', function() {
                const isValid = validateSequence(this.value);
                
                validator.className = `w-3 h-3 rounded-full ${isValid ? 'bg-green-500' : 'bg-red-500'}`;
                text.textContent = isValid ? 'Valid' : 'Invalid';
                text.className = isValid ? 'text-green-600 semi-bold dark:text-green-400' : 'text-red-600 dark:text-red-400';
            });
        });
    </script>
</body>
</html>
    ''', result=result)         

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=False)