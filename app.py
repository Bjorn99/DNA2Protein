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
app.config['API_KEY'] = os.getenv('API_KEY', '83469fae42ffbc3846e3adbc8bc02fb09609')

# NCBI access
Entrez.api_key = app.config('API_KEY')

class DNAValidator:
    """Validates DNA sequences and provides basic DNA operations."""
    
    @staticmethod
    def validate_sequence(sequence: str) -> bool:
        """
        Validate if the sequence contains only valid DNA nucleotides.
        """
        valid_nucleotides = set('ATCG')
        return all(nucleotide in valid_nucleotides for nucleotide in sequence.upper())

    @staticmethod
    def clean_sequence(sequence: str) -> str:
        """
        Clean and standardize DNA sequence.
        """
        return ''.join(sequence.split()).upper()

    @staticmethod
    def check_sequence_quality(sequence: str) -> tuple[bool, str]:
        """
        Check sequence quality and return validation status and message.
        """
        if not sequence:
            return False, "Empty sequence provided"
        if len(sequence) < 3:
            return False, "Sequence too short (minimum 3 nucleotides required)"
        if len(set(sequence)) == 1:
            return False, "Sequence consists of a single repeated nucleotide"
        return True, "Sequence valid"

class GeneticCode:
    """Manages genetic code and codon frequencies."""
    
    def __init__(self):
        self.code = self._initialize_genetic_code()
        
    def _initialize_genetic_code(self) -> Dict:
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

        return gencode
     def get_amino_acid(self, codon: str) -> str:
        """Get amino acid for a given codon."""
        return self.code.get(codon, {}).get('aa', 'X')
    
    def get_frequency(self, codon: str) -> float:
        """Get usage frequency for a given codon."""
        return self.code.get(codon, {}).get('freq', 0.0)

class SequenceAnalyzer:
    """Analyzes DNA sequences for various features."""
    
    def __init__(self):
        self.genetic_code = GeneticCode()
        self.validator = DNAValidator()
        self.kozak_regex = re.compile(r'(G|A)NN(A|G)TGATG')
        
    def find_orfs(self, sequence: str) -> List[str]:
        """Find all open reading frames in the sequence."""
        pattern = re.compile(r'(?=(ATG(?:...)*?(?:TAA|TAG|TGA)))')
        return pattern.findall(sequence)
    
    def translate_dna(self, sequence: str) -> str:
        """Translate DNA sequence to protein."""
        protein = []
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon in {'TAA', 'TAG', 'TGA'}:
                break
            amino_acid = self.genetic_code.get_amino_acid(codon)
            protein.append(amino_acid)
        return ''.join(protein)

    def find_kozak_sequences(self, sequence: str) -> List[int]:
        """Find Kozak consensus sequences."""
        return [match.start() for match in self.kozak_regex.finditer(sequence)]

    def calculate_cai(self, sequence: str) -> float:
        """Calculate Codon Adaptation Index."""
        codons = [sequence[i:i+3] for i in range(0, len(sequence) - 2, 3)]
        if not codons:
            return 0.0
        
        frequencies = [self.genetic_code.get_frequency(codon) for codon in codons if codon in self.genetic_code.code]
        if not frequencies:
            return 0.0
            
        return sum(frequencies) / len(frequencies)


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

    def predict(self, sequence: str) -> Dict:
        """Predict signal peptide presence and characteristics."""
        if len(sequence) < 30:
            return {"error": "Sequence too short for signal peptide prediction"}
            
        features = self._extract_features(sequence[:30])
        prediction = self.model.predict_proba(features)[0]
        
        return {
            'is_signal_peptide': bool(prediction.argmax()),
            'confidence': float(max(prediction) * 100),
            'details': {
                'hydrophobicity_profile': self._calculate_hydrophobicity_profile(sequence[:30]),
                'charge_distribution': self._calculate_charge_distribution(sequence[:30])
            }
        }

    def _extract_features(self, sequence: str) -> np.ndarray:
        """Extract sequence features for prediction."""
        features = []
        window_size = 20

        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]
            hydrophobicity = sum(self.amino_acid_properties[aa]['hydrophobicity'] for aa in window)
            charge = sum(self.amino_acid_properties[aa]['charge'] for aa in window)
            size = sum(self.amino_acid_properties[aa]['size'] for aa in window)
            features.extend([hydrophobicity, charge, size])
            
        return np.array(features).reshape(1, -1)

def analyze_dna(dna_sequence: str, include_extended: bool = False) -> Dict:
    """Main function to analyze DNA sequence."""
    
    analyzer = SequenceAnalyzer()
    validator = DNAValidator()
    
    # Clean and validate sequence
    sequence = validator.clean_sequence(dna_sequence)
    is_valid, message = validator.check_sequence_quality(sequence)
    
    if not is_valid:
        return {"error": message}
    
    if not validator.validate_sequence(sequence):
        return {"error": "Invalid DNA sequence. Please use only A, T, C, and G."}
        
    try:
        # Find ORFs
        orfs = analyzer.find_orfs(sequence)
        if not orfs:
            return {"error": "No open reading frames found."}
            
        # Get longest ORF and translate
        longest_orf = max(orfs, key=len)
        protein = analyzer.translate_dna(longest_orf)
        
        # Basic analysis
        result = {
            "longest_orf": longest_orf,
            "protein": protein,
            "kozak_positions": analyzer.find_kozak_sequences(sequence),
            "cai": analyzer.calculate_cai(longest_orf),
            "gc_content": GC(sequence),
            "molecular_weight": molecular_weight(sequence),
            "protein_length": len(protein)
        }
        
        # Signal peptide analysis
        signal_peptide_predictor = SignalPeptidePredictor()
        result["signal_peptide"] = signal_peptide_predictor.predict(protein)
        
        # Extended analysis if requested
        if include_extended:
            db_integrator = BiologicalDatabaseIntegrator()
            result.update({
                "similar_sequences": db_integrator.fetch_similar_sequences(protein),
                "protein_domains": db_integrator.get_protein_domains(protein)
            })
            
        return result
        
    except Exception as e:
        return {"error": f"Analysis error: {str(e)}"}


@app.route('/', methods=['GET', 'POST'])
def index():
    result = None
    if request.method == "POST":
        dna_sequence = request.form['dna_sequence']
        include_extended = request.form.get('include_extended', False)
        result = analyze_dna(dna_sequence, include_extended)

    return render_template_string('''   
        <!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>DNA2Protein - DNA Analysis Tool</title>
    <link href="https://cdnjs.cloudflare.com/ajax/libs/tailwindcss/2.2.19/tailwind.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0-beta3/css/all.min.css" rel="stylesheet">
    <script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.24.2/plotly.min.js"></script>
    
    <style>
        /* Base Styles */
        body {
            background: linear-gradient(135deg, #f0f4f8 0%, #e2e8f0 100%);
            color: #2d3748;
            font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
            min-height: 100vh;
        }

        /* Container */
        .container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 2rem;
        }

        /* Title */
        .title {
            color: #2b6cb0;
            font-size: 3rem;
            font-weight: 800;
            text-align: center;
            margin-bottom: 2rem;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.1);
            letter-spacing: -0.025em;
        }

        /* Feature Grid */
        .feature-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 1.5rem;
            margin-bottom: 2rem;
        }

        .feature-card {
            background: white;
            border-radius: 1rem;
            padding: 1.5rem;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
            transition: all 0.3s ease;
        }

        .feature-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 10px 15px rgba(0, 0, 0, 0.1);
        }

        /* Input Form */
        .input-form {
            background: white;
            padding: 2rem;
            border-radius: 1rem;
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);
            margin-bottom: 2rem;
        }

        .input-label {
            display: block;
            font-weight: 600;
            margin-bottom: 0.5rem;
            color: #4a5568;
        }

        .input-field {
            width: 100%;
            padding: 0.75rem;
            border: 2px solid #e2e8f0;
            border-radius: 0.5rem;
            font-size: 1rem;
            transition: all 0.3s ease;
            margin-bottom: 1.5rem;
        }

        .input-field:focus {
            border-color: #4299e1;
            box-shadow: 0 0 0 3px rgba(66, 153, 225, 0.5);
            outline: none;
        }

        /* Submit Button */
        .submit-button {
            background-color: #4299e1;
            color: white;
            font-weight: 600;
            padding: 0.75rem 2rem;
            border-radius: 0.5rem;
            transition: all 0.3s ease;
            display: flex;
            align-items: center;
            justify-content: center;
            gap: 0.5rem;
        }

        .submit-button:hover {
            background-color: #3182ce;
            transform: translateY(-2px);
        }

        .submit-button:active {
            transform: translateY(0);
        }

        /* Results Section */
        .result-box {
            background: white;
            border-radius: 1rem;
            padding: 2rem;
            margin-top: 2rem;
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);
        }

        .result-title {
            font-size: 1.5rem;
            font-weight: 700;
            color: #2b6cb0;
            margin-bottom: 1.5rem;
            padding-bottom: 0.5rem;
            border-bottom: 2px solid #e2e8f0;
        }

        .result-item {
            margin-bottom: 1.5rem;
        }

        .result-label {
            font-weight: 600;
            color: #4a5568;
            margin-bottom: 0.5rem;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }

        .result-value {
            background: #f7fafc;
            padding: 1rem;
            border-radius: 0.5rem;
            font-family: 'Courier New', monospace;
            word-break: break-all;
            line-height: 1.5;
        }

        /* Loading Overlay */
        .loading-overlay {
            position: fixed;
            inset: 0;
            background: rgba(255, 255, 255, 0.9);
            display: flex;
            justify-content: center;
            align-items: center;
            z-index: 50;
        }

        .loading-spinner {
            width: 50px;
            height: 50px;
            border: 4px solid #e2e8f0;
            border-top-color: #4299e1;
            border-radius: 50%;
            animation: spin 1s linear infinite;
        }

        @keyframes spin {
            to { transform: rotate(360deg); }
        }

        /* Error Message */
        .error-message {
            background: #fff5f5;
            color: #c53030;
            padding: 1rem;
            border-radius: 0.5rem;
            border: 1px solid #feb2b2;
            margin-top: 1rem;
        }

        /* Footer */
        .footer {
            text-align: center;
            margin-top: 2rem;
            padding: 1rem;
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

        /* Responsive Adjustments */
        @media (max-width: 640px) {
            .container {
                padding: 1rem;
            }

            .title {
                font-size: 2rem;
            }

            .feature-grid {
                grid-template-columns: 1fr;
            }

            .result-value {
                font-size: 0.875rem;
            }
        }

        /* Copy Button */
        .copy-button {
            padding: 0.25rem 0.5rem;
            border-radius: 0.25rem;
            font-size: 0.875rem;
            color: #4299e1;
            transition: all 0.2s ease;
        }

        .copy-button:hover {
            background-color: #ebf8ff;
            color: #2b6cb0;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1 class="title">DNA2Protein</h1>
        
        <!-- Feature Grid -->
        <div class="feature-grid">
            <div class="feature-card">
                <h3 class="text-xl font-bold text-blue-600 mb-2">ORF Analysis</h3>
                <p class="text-gray-600">Identifies open reading frames in DNA sequences and predicts protein coding regions</p>
            </div>
            <div class="feature-card">
                <h3 class="text-xl font-bold text-blue-600 mb-2">Protein Translation</h3>
                <p class="text-gray-600">Translates DNA sequences to amino acid sequences using standard genetic code</p>
            </div>
            <div class="feature-card">
                <h3 class="text-xl font-bold text-blue-600 mb-2">Signal Peptide Analysis</h3>
                <p class="text-gray-600">Predicts presence and characteristics of signal peptides in protein sequences</p>
            </div>
        </div>

        <!-- Input Form -->
        <form method="post" class="input-form">
            <label for="dna_sequence" class="input-label">Enter DNA Sequence:</label>
            <textarea id="dna_sequence" name="dna_sequence" required class="input-field" 
                      placeholder="e.g., ATGCGATCGATCG..." rows="4"></textarea>
            <div class="flex items-center gap-4">
                <button type="submit" class="submit-button">
                    <i class="fas fa-dna"></i>
                    Analyze DNA
                </button>
                <label class="flex items-center gap-2">
                    <input type="checkbox" name="include_extended" class="form-checkbox h-4 w-4 text-blue-600">
                    <span class="text-sm text-gray-600">Include extended analysis</span>
                </label>
            </div>
        </form>

        <!-- Loading Overlay -->
        <div id="loading-overlay" class="loading-overlay hidden">
            <div class="loading-spinner"></div>
        </div>

        <!-- Results Section -->
        {% if result %}
        <div class="result-box">
            <h2 class="result-title">Analysis Results</h2>
            {% if result.error %}
            <div class="error-message">{{ result.error }}</div>
            {% else %}
            <div class="result-item">
                <div class="result-label">
                    <span>Longest Open Reading Frame (ORF)</span>
                    <button onclick="copyToClipboard('orf-sequence')" class="copy-button">
                        <i class="fas fa-copy"></i> Copy
                    </button>
                </div>
                <div id="orf-sequence" class="result-value">{{ result.longest_orf }}</div>
            </div>
            <div class="result-item">
                <div class="result-label">
                    <span>Translated Protein</span>
                    <button onclick="copyToClipboard('protein-sequence')" class="copy-button">
                        <i class="fas fa-copy"></i> Copy
                    </button>
                </div>
                <div id="protein-sequence" class="result-value">{{ result.protein }}</div>
            </div>
            <div class="result-item">
                <div class="result-label">Codon Adaptation Index (CAI)</div>
                <div class="result-value">{{ "%.3f"|format(result.cai) }}</div>
            </div>
            <div class="result-item">
                <div class="result-label">GC Content</div>
                <div class="result-value">{{ "%.1f"|format(result.gc_content) }}%</div>
            </div>
            <div class="result-item">
                <div class="result-label">Molecular Weight</div>
                <div class="result-value">{{ "%.2f"|format(result.molecular_weight) }} Da</div>
            </div>
            <div class="result-item">
                <div class="result-label">Signal Peptide Prediction</div>
                <div class="result-value">
                    <div>Probability: {{ "%.1f"|format(result.signal_peptide.confidence) }}%</div>
                    <div>Prediction: {{ "Signal peptide detected" if result.signal_peptide.is_signal_peptide else "No signal peptide detected" }}</div>
                </div>
            </div>
            {% endif %}
        </div>
        {% endif %}

        <!-- Footer -->
        <div class="footer">
            <p>Created by <a href="https://github.com/Bjorn99" target="_blank">Bjorn99</a></p>
            <p>Have issues or want to contribute? Visit the <a href="https://github.com/Bjorn99/DNA2Protein" target="_blank">GitHub repository</a></p>
        </div>
    </div>

    <script>
        // Copy to clipboard functionality
        function copyToClipboard(elementId) {
            const element = document.getElementById(elementId);
            const text = element.textContent;
            navigator.clipboard.writeText(text).then(() => {
                const button = element.previousElementSibling.querySelector('.copy-button');
                const originalHTML = button.innerHTML;
                button.innerHTML = '<i class="fas fa-check"></i> Copied!';
                setTimeout(() => {
                    button.innerHTML = originalHTML;
                }, 2000);
            });
        }

        // Form submission handling
        document.querySelector('form').addEventListener('submit', function() {
            document.getElementById('loading-overlay').classList.remove('hidden');
        });

        // Initialize tooltips if using them
        if (typeof tippy !== 'undefined') {
            tippy('[data-tippy-content]');
        }

        // Optional: Add sequence validation
        document.getElementById('dna_sequence').addEventListener('input', function(e) {
            const sequence = e.target.value.toUpperCase();
            const validNucleotides = new Set(['A', 'T', 'C', 'G']);
            const isValid = [...sequence].every(char => 
                validNucleotides.has(char) || char === ' ' || char === '\n'
            );
            
            e.target.classList.toggle('border-red-500', !isValid);
            const submitButton = document.querySelector('.submit-button');
            submitButton.disabled = !isValid;
            submitButton.classList.toggle('opacity-50', !isValid);
        });
    </script>
</body>
</html>
    ''', result=result)         

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=True)