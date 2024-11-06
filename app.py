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

def calculate_cai(dna):
    """Calculate the Codon Adaptation Index (CAI) for a DNA sequence."""
    codon_counts = Counter(dna[i:i+3] for i in range(0, len(dna) - 2, 3))
    total_codons = sum(codon_counts.values())
    return sum(count / total_codons * len(gencode[codon]) for codon, count in codon_counts.items() if codon in gencode) / (len(dna) // 3)

def predict_signal_peptide(protein):
    """Simple prediction of signal peptide presence based on N-terminal sequence."""
    n_terminal = protein[:30]  # Consider first 30 amino acids
    if n_terminal.count('L') + n_terminal.count('A') + n_terminal.count('V') > 10:
        return "Potential signal peptide detected"
    return "No signal peptide detected"

def analyze_dna(dna):
    """Analyze DNA sequence for ORFs, Kozak sequences, and translate to protein."""
    dna = ''.join(dna.split()).upper()
    
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
    
    return {
        "longest_orf": longest_orf,
        "protein": protein,
        "kozak_positions": kozak_positions,
        "cai": cai,
        "signal_peptide": signal_peptide
    }

@app.route('/', methods=['GET', 'POST'])
def index():
    result = None
    if request.method == "POST":
        dna_sequence = request.form['dna_sequence']
        result = analyze_dna(dna_sequence)

    return render_template_string('''   
        <!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>DNA to Protein Translator</title>
    <link href="https://cdn.jsdelivr.net/npm/tailwindcss@2.2.19/dist/tailwind.min.css" rel="stylesheet">
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0-beta3/css/all.min.css">
    <style>
        /* Base Styles */
        body {
            background: linear-gradient(135deg, #f6f8fa 0%, #e9ecef 100%);
            color: #2d3748;
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            min-height: 100vh;
        }

        /* Container */
        .container {
            max-width: 800px;
            margin: 0 auto;
            padding: 2rem;
        }

        /* Title Styles */
        .title {
            color: #2C3E50;
            font-size: 2.75rem;
            font-weight: 800;
            text-align: center;
            margin-bottom: 2rem;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.1);
            position: relative;
            padding-bottom: 1rem;
        }

        .title::after {
            content: '';
            position: absolute;
            bottom: 0;
            left: 50%;
            transform: translateX(-50%);
            width: 100px;
            height: 4px;
            background: linear-gradient(90deg, #2C3E50 0%, #34495E 100%);
            border-radius: 2px;
        }

        /* Form Styles */
        .input-form {
            background-color: white;
            padding: 2.5rem;
            border-radius: 16px;
            box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.1);
            backdrop-filter: blur(10px);
            border: 1px solid rgba(255, 255, 255, 0.2);
        }

        .input-label {
            font-weight: 600;
            margin-bottom: 0.75rem;
            color: #2C3E50;
            font-size: 1.2rem;
            display: block;
        }

        .input-field {
            width: 100%;
            padding: 1rem;
            border: 2px solid #e2e8f0;
            border-radius: 12px;
            font-size: 1.1rem;
            transition: all 0.3s ease;
            background-color: #f8fafc;
            margin-bottom: 1.5rem;
        }

        .input-field:focus {
            border-color: #2C3E50;
            outline: none;
            box-shadow: 0 0 0 3px rgba(44, 62, 80, 0.2);
            background-color: white;
        }

        /* Button Styles */
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

        /* Result Box Styles */
        .result-box {
            background-color: white;
            border-radius: 16px;
            padding: 2rem;
            margin-top: 2rem;
            box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.1);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: all 0.3s ease;
        }

        .result-box:hover {
            transform: translateY(-5px);
            box-shadow: 0 20px 30px -10px rgba(0, 0, 0, 0.15);
        }

        .result-title {
            font-size: 1.75rem;
            font-weight: 700;
            color: #2C3E50;
            margin-bottom: 1.5rem;
            padding-bottom: 1rem;
            border-bottom: 2px solid #e2e8f0;
            position: relative;
        }

        .result-title::after {
            content: '';
            position: absolute;
            bottom: -2px;
            left: 0;
            width: 50px;
            height: 2px;
            background: #2C3E50;
        }

        .result-item {
            margin-bottom: 1.5rem;
            transition: all 0.3s ease;
            padding: 1rem;
            border-radius: 12px;
            background-color: #f8fafc;
        }

        .result-item:hover {
            transform: translateX(5px);
            background-color: #f1f5f9;
        }

        .result-label {
            font-weight: 600;
            color: #2C3E50;
            margin-bottom: 0.5rem;
            font-size: 1.1rem;
        }

        .result-value {
            background-color: white;
            padding: 1rem;
            border-radius: 8px;
            font-family: 'Courier New', Courier, monospace;
            word-break: break-all;
            font-size: 1rem;
            border: 1px solid #e2e8f0;
            transition: all 0.3s ease;
        }

        .result-value:hover {
            border-color: #2C3E50;
            box-shadow: 0 2px 8px rgba(44, 62, 80, 0.1);
        }

        /* Error Message Styles */
        .error-message {
            color: #dc2626;
            font-weight: 600;
            margin-top: 1rem;
            padding: 1rem;
            border-radius: 12px;
            background-color: #fef2f2;
            border: 1px solid #fee2e2;
            animation: shake 0.5s ease-in-out;
        }

        /* Footer Styles */
        .footer {
            margin-top: 3rem;
            text-align: center;
            font-size: 1rem;
            color: #4b5563;
            padding: 2rem 0;
            border-top: 1px solid #e5e7eb;
        }

        .footer a {
            color: #2C3E50;
            text-decoration: none;
            font-weight: 600;
            transition: color 0.3s ease;
        }

        .footer a:hover {
            color: #34495E;
            text-decoration: underline;
        }

        /* Animations */
        @keyframes dna-spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }

        @keyframes shake {
            0%, 100% { transform: translateX(0); }
            25% { transform: translateX(-5px); }
            75% { transform: translateX(5px); }
        }

        /* Responsive Design */
        @media (max-width: 640px) {
            .container {
                padding: 1rem;
            }

            .title {
                font-size: 2rem;
            }

            .input-form {
                padding: 1.5rem;
            }

            .result-box {
                padding: 1.5rem;
            }
        }

        /* Dark Mode Support */
        @media (prefers-color-scheme: dark) {
            body {
                background: linear-gradient(135deg, #1a202c 0%, #2d3748 100%);
                color: #e2e8f0;
            }

            .input-form, .result-box {
                background-color: #2d3748;
                border-color: rgba(255, 255, 255, 0.1);
            }

            .input-field {
                background-color: #1a202c;
                border-color: #4a5568;
                color: #e2e8f0;
            }

            .result-item {
                background-color: #1a202c;
            }

            .result-value {
                background-color: #2d3748;
                border-color: #4a5568;
                color: #e2e8f0;
            }

            .title, .result-title, .result-label {
                color: #e2e8f0;
            }

            .footer {
                border-color: #4a5568;
                color: #e2e8f0;
            }

            .footer a {
                color: #90cdf4;
            }
        }
    </style>
</head>
<body class="min-h-screen">
    <div class="container">
        <h1 class="title">DNA2PROTEIN</h1>
        <form method="post" class="input-form">
            <label for="dna_sequence" class="input-label">Enter DNA Sequence:</label>
            <input type="text" 
                   id="dna_sequence" 
                   name="dna_sequence" 
                   required 
                   class="input-field" 
                   placeholder="e.g., ATGCGATCGATCG">
            <button type="submit" class="submit-button">
                <i class="fas fa-dna"></i> Analyze DNA
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
        </div>
        {% endif %}

        <div class="footer">
            <p class="mb-2">Created by <a href="https://github.com/Bjorn99" target="_blank">Bjorn99</a></p>
            <p>Have issues or want to contribute? Visit the <a href="https://github.com/Bjorn99/DNA2Protein" target="_blank">GitHub repository</a></p>
        </div>
    </div>
</body>
</html>
    ''', result=result)         

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=True)