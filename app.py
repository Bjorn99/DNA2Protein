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
        "reverse_complement": reverse_complement(dna)
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
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>DNA2Protein Analysis Tool</title>
    <link href="https://cdn.jsdelivr.net/npm/tailwindcss@2.2.19/dist/tailwind.min.css" rel="stylesheet">
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0-beta3/css/all.min.css">
    <style>
        :root {
            --primary-color: #2C3E50;
            --secondary-color: #34495E;
            --accent-color: #576574;
            --success-color: #27ae60;
            --error-color: #e74c3c;
            --text-primary: #1a202c;
            --text-secondary: #2d3748;
            --text-light: #718096;
            --text-dark: #2d3748;
        }

        body {
            background: linear-gradient(135deg, #f6f8fa 0%, #e9ecef 100%);
            color: #2d3748;
            font-family: 'Inter', system-ui, -apple-system, sans-serif;
        }

        .container {
            max-width: 1000px;
            margin: 2rem auto;
            padding: 0 1.5rem;
        }

        .title {
            font-size: 3.5rem;
            font-weight: 800;
            background: linear-gradient(135deg, #1a202c 0%, #2d3748 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            text-align: center;
            margin: 2rem 0;
            letter-spacing: -0.05em;
        }

        .input-form {
            background: #2d3748;
            bordere: 1px solid #4a5568;
            padding: 2rem;
            border-radius: 1rem;
            box-shadow: 0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04);
            backdrop-filter: blur(10px);
            margin-bottom: 2rem;
        }

        .input-field {
            background: #1a202c;
            border-color: #4a5568;
            color: var(--text-primary);
            width: 100%;
            padding: 1rem;
            border: 2px solid #E5E7EB;
            border-radius: 0.75rem;
            font-size: 1.1rem;
            transition: all 0.2s ease;
            margin-bottom: 1.5rem;
            font-family: 'Courier New', monospace;
        }

        .input-field:focus {
            border-color: var(--primary-color);
            box-shadow: 0 0 0 3px rgba(44, 62, 80, 0.2);
            outline: none;
        }

        .submit-button {
            position: relative;
            overflow: hidden;
            background: linear-gradient(135deg, #2C3E50 0%, #34495E 100%);
            color: white;
            font-weight: 600;
            padding: 1rem 2.5rem;
            border-radius: 12px;
            font-size: 1.1rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            transition: all 0.3s ease;
            width: 100%;
            border: none;
            cursor: pointer;
        }

        .submit-button:hover {
            transform: translateY(-2px);
            box-shadow: 0 8px 15px rgba(44, 62, 80, 0.3);
            background: linear-gradient(135deg, #34495E 0%, #2C3E50 100%);
        }

        .submit-button i {
            margin-right: 0.75rem;
            transition: transform 0.3s ease;
        }

        .submit-button:hover i {
            animation: dna-spin 2s linear infinite;
        }

        @keyframes dna-spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }

        .result-box {
            background: #2d3748;
            border: 1px solid #4a5568;
            border-radius: 1rem;
            padding: 2rem;
            margin-top: 2rem;
            box-shadow: 0 20px 25px -5px rgba(0, 0, 0, 0.1);
            transition: transform 0.3s ease;
        }

        .result-box:hover {
            transform: translateY(-5px);
        }

        .result-item {
            background: #1a202c;
            border-radius: 0.75rem;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
            transition: all 0.2s ease;
        }

        .result-item:hover {
            background: #F9FAFB;
            transform: scale(1.01);
        }

        .result-label {
            backgroung: #2d3748;
            color: #1a202c;
            font-weight: 600;
            margin-bottom: 0.5rem;
            font-size: 1.1rem;
        }

        .result-label:hover {
            color: black;
        }

        .result-value {
            font-family: 'Courier New', monospace;
            background: #4a5568;
            color: var(--text-primary);
            padding: 1rem;
            border-radius: 0.5rem;
            word-break: break-all;
            border: 1px solid #E5E7EB;
        }

        .sequence-validator {
            color: var(--text-secondary);s
            display: flex;
            align-items: center;
            gap: 0.5rem;
            margin-top: 0.5rem;
            margin-bottom: 1rem;
            font-size: 0.9rem;
        }

        .validator-icon {
            width: 16px;
            height: 16px;
            border-radius: 50%;
        }

        .validator-icon.valid {
            background-color: var(--success-color);
        }

        .validator-icon.invalid {
            background-color: var(--error-color);
        }

        /* Responsive Design */
        @media (max-width: 768px) {
            .container {
                margin: 1rem auto;
            }

            .title {
                font-size: 2.5rem;
            }

            .input-form, .result-box {
                padding: 1.5rem;
            }
        }

        /* Dark Mode Improvements */
        @media (prefers-color-scheme: dark) {
            :root {
                --text-primary: #f7fafc;
                --text-secondary: #e2e8f0;
                --text-light: #cbd5e0;
            }

            body {
                background: linear-gradient(135deg, #1a202c 0%, #2d3748 100%);
                color: var(--text-primary);
            }

            .title {
                /* Lighter gradient for dark mode */
                background: linear-gradient(135deg, #f7fafc 0%, #e2e8f0 100%);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
            }

            .input-form, .result-box {
                background: #2d3748;
                border: 1px solid #4a5568;
            }

            .input-field {
                background: #1a202c;
                border-color: #4a5568;
                color: var(--text-primary);
            }

            .result-item {
                background: #1a202c;
            }

            .result-label {
                color: #e2e8f0;  /* Lighter color for dark mode */
            }

            .result-value {
                background: #2d3748;
                border-color: #4a5568;
                color: var(--text-primary);
            }

            /* Footer links in dark mode */
            footer a {
                color: #90cdf4 !important;  /* Lighter blue for better visibility */
            }

            footer a:hover {
                color: green !important;
            }

            /* Validator text colors for dark mode */
            .sequence-validator {
                color: var(--text-light);
            }

            .validator-text.valid {
                color: #68d391 !important;  /* Lighter green for dark mode */
            }

            .validator-text.invalid {
                color: #fc8181 !important;  /* Lighter red for dark mode */
            }
        }

        /* Error message improvements */
        .bg-red-100 {
            background-color: #fff5f5;
            color: #c53030;
        }

        @media (prefers-color-scheme: dark) {
            .bg-red-100 {
                background-color: #742a2a;
                color: #feb2b2;
                border-color: #fc8181;
            }
        }

    </style>
</head>
<body>
    <div class="container">
        <h1 class="title">DNA2PROTEIN</h1>
        
        <div class="input-form">
            <form method="post" id="dna-form">
                <label for="dna_sequence" class="block text-lg font-semibold mb-2">Enter DNA Sequence:</label>
                <input type="text" 
                       id="dna_sequence" 
                       name="dna_sequence" 
                       required 
                       class="input-field" 
                       placeholder="e.g., ATGCGATCGATCG"
                       title="Please enter valid DNA sequence (A, T, C, G only)">
                
                <div class="sequence-validator">
                    <div class="validator-icon"></div>
                    <span class="validator-text">Enter a valid DNA sequence</span>
                </div>

                <button type="submit" class="submit-button">
                    <i class="fas fa-dna"></i>
                    Analyze DNA
                </button>
            </form>
        </div>

        {% if error %}
        <div class="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded relative" role="alert">
            <strong class="font-bold">Error: </strong>
            <span class="block sm:inline">{{ error }}</span>
        </div>
        {% endif %}

        {% if result %}
        <div class="result-box">
            <h2 class="text-2xl font-bold mb-6 pb-2 border-b-2 border-gray-200">Analysis Results</h2>
            
            {% if result.error %}
            <div class="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded relative" role="alert">
                <strong class="font-bold">Error: </strong>
                <span class="block sm:inline">{{ result.error }}</span>
            </div>
            {% else %}
            <div class="result-item">
                <div class="result-label">
                    <i class="fas fa-dna mr-2"></i>Longest Open Reading Frame (ORF):
                </div>
                <div class="result-value">{{ result.longest_orf }}</div>
            </div>
            
            <div class="result-item">
                <div class="result-label">
                    <i class="fas fa-project-diagram mr-2"></i>Translated Protein:
                </div>
                <div class="result-value">{{ result.protein }}</div>
            </div>
            
            <div class="result-item">
                <div class="result-label">
                    <i class="fas fa-map-marker-alt mr-2"></i>Kozak Sequence Positions:
                </div>
                <div class="result-value">
                    {% if result.kozak_positions %}
                        {{ result.kozak_positions|join(', ') }}
                    {% else %}
                        No Kozak sequences found
                    {% endif %}
                </div>
            </div>
            
            <div class="result-item">
                <div class="result-label">
                    <i class="fas fa-chart-line mr-2"></i>Codon Adaptation Index (CAI):
                </div>
                <div class="result-value">{{ "%.2f"|format(result.cai) }}</div>
            </div>
            
            <div class="result-item">
                <div class="result-label">
                    <i class="fas fa-microscope mr-2"></i>Signal Peptide Prediction:
                </div>
                <div class="result-value">{{ result.signal_peptide }}</div>
            </div>

            <div class="result-item">
                <div class="result-label">
                    <i class="fas fa-percentage mr-2"></i>GC Content:
                </div>
                <div class="result-value">{{ result.gc_content }}%</div>
            </div>

            <div class="result-item">
                <div class="result-label">
                    <i class="fas fa-calculator mr-2"></i>Nucleotide Frequency:
                </div>
                <div class="result-value">
                    A: {{ result.nucleotide_frequency.A }},
                    T: {{ result.nucleotide_frequency.T }},
                    G: {{ result.nucleotide_frequency.G }},
                    C: {{ result.nucleotide_frequency.C }}
                </div>
            </div>

            <div class="result-item">
                <div class="result-label">
                    <i class="fas fa-ruler mr-2"></i>Sequence Length:
                </div>
                <div class="result-value">{{ result.sequence_length }} bp</div>
            </div>

            <div class="result-item">
                <div class="result-label">
                    <i class="fas fa-arrow-left mr-2"></i>Reverse Complement:
                </div>
                <div class="result-value">{{ result.reverse_complement }}</div>
            </div>
            {% endif %}
        </div>
        {% endif %}

        <footer class="mt-8 text-center">
            <p class="mb-2">Created by <a href="https://github.com/Bjorn99" class="hover:text-green-900 dark:text-white-400" target="_blank">Bjorn99</a></p>
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
            const form = document.getElementById('dna-form');
            const input = document.getElementById('dna_sequence');
            const validator = document.querySelector('.sequence-validator');
            const validatorIcon = document.querySelector('.validator-icon');
            const validatorText = document.querySelector('.validator-text');

            function validateSequence(sequence) {
                const validChars = /^[ATCGatcg\\s]+$/;
                return validChars.test(sequence);
            }

            input.addEventListener('input', function() {
                const sequence = this.value;
                const isValid = validateSequence(sequence);
                
                validatorIcon.className = 'validator-icon ' + (isValid ? 'valid' : 'invalid');
                validatorText.textContent = isValid ? 'Valid DNA sequence' : 'Invalid characters detected';
                validatorText.style.color = isValid ? 'var(--success-color)' : 'var(--error-color)';
            });
        });
    </script>
</body>
</html>
    ''', result=result)         

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=False)