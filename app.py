from flask import Flask, request, render_template_string
from dotenv import load_dotenv
import os

# Load environment variables from .env file
load_dotenv()

# Create a Flask application instance
app = Flask(__name__)

# Configure the application with environment variables
app.config['SECRET_KEY'] = os.getenv('SECRET_KEY')  # Secret key for sessions and CSRF protection
app.config['DATABASE_URL'] = os.getenv('DATABASE_URL')  # Database connection URL (if needed)
app.config['API_KEY'] = os.getenv('API_KEY')  # API key for external services (if needed)


app = Flask(__name__)


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
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}

# STOP CODONS "TAA", "TAG", "TGA"

def translate_dna(dna):
    # Remove spaces from the DNA sequence
    dna = dna.replace(" ", "").replace("\n", "").upper( )

    # check if the length of dna is a multiple of 3 
    if len(dna) % 3 != 0:
        dna = dna[:len(dna) - (len(dna) % 3)] # Truncate to nearest multiple of 3
        
    protein = []
    stop_codons = {'TAA', 'TAG', 'TGA'} 

    for start in range(0, len(dna), 3):
        codon = dna[start:start+3]
        aa = gencode.get(codon.upper(), "X")  # Use "X" for unknown codons
        protein.append(aa)

        if codon in stop_codons:
            break
    return ''.join(protein)

@app.route('/', methods=['GET', 'POST'])
def index():
    protein = ""
    if request.method == "POST":
        dna_sequence = request.form['dna_sequence']
        protein = translate_dna(dna_sequence)
    
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
             body {
            background-color: #f7fafc; /* Light grayish-blue background */
        }

        /* Centering the container */
        .container {
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            min-height: 100vh; /* Full height */
        }

        input[type="text"] {
            width: auto; /* Make width auto to flex according to content */
            max-width: 100%; /* Prevent overflow beyond container */
            padding: 10px; /* Add padding for better spacing */
            border: 2px solid #38a169; /* Green border */
            border-radius: 5px; /* Rounded corners */
            transition: border-color 0.3s ease; /* Smooth transition for focus */
        }

        input[type="text"]:focus {
            border-color: #2f855a; /* Darker green on focus */
            outline: none; /* Remove default outline */
        }   

        /* Styling for the result box */
        .result-box {
    background-color: #f0fff4; /* Light green background */
    border: 2px solid #38a169; /* Green border */
    border-radius: 8px; /* Rounded corners */
    padding: 20px; /* Padding around text */
    width: auto; /* Width adjusts based on content */
    max-width: 90%; /* Maximum width */
    transition: all 0.3s ease; /* Smooth transition */
    opacity: 0; /* Start hidden for animation */
    transform: translateY(-20px); /* Slide up effect */
    white-space: pre-wrap; /* Preserve whitespace and wrap text */
    word-wrap: break-word; /* Break long words onto the next line */
    overflow-wrap: break-word; /* Break words if they are too long */
}

        /* Animation when showing result */
        .result-box.show {
            opacity: 1; /* Fully visible */
            transform: translateY(0); /* Move to original position */
        }

        #arrow-icon {
            transition: transform 0.3s ease; /* Smooth transition for movement */
        }

        #submit-btn:hover #arrow-icon {
            transform: translateX(5px); /* Move the arrow to the right on hover */
        }
    </style>



        </head>
        <body class="bg-gray-100">
            <div class="container mx-auto p-4">
                <h1 class="text-3xl font-bold text-center mb-4">Translate DNA to Protein</h1>
                <form method="post" class="bg-white p-6 rounded shadow-md">
                    <label for="dna_sequence" class="block text-gray-700 flex">Enter DNA Sequence:</label>
                    <input type="text" id="dna_sequence" name="dna_sequence" required class="border border-gray-300 p-2 w-full rounded mb-4">
                    
                     <!-- Centering the button -->
                    <div class="flex justify-center mb-4">
                        <button id="submit-btn" type="submit" class="bg-blue-500 hover:bg-green-600 text-white font-bold py-1.5 px-10 rounded-md flex items-center transition-transform duration-300"><i class="fas fa-arrow-right mr-2 transform transition-transform duration-300 hover:translate-x-1" id="arrow-icon"></i>
                        Translate
                        </button>
                    </div>
                </form>
                
                {% if protein %}
                    <h2 class="text-xl font-semibold mt-4 text-center">Protein Sequence:</h2>
                    <div class="result-box show">{{ protein }}</div>
                {% endif %}
            </div>

            <script>
        // Optional JavaScript to add animation class after rendering
        document.addEventListener("DOMContentLoaded", function() {
            const resultBox = document.querySelector('.result-box');
            if (resultBox) {
                resultBox.classList.add('show'); // Adding show class to trigger animation
            }
        });
    </script>


        </body>
        </html>
    ''', protein=protein)

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))  # Get the port from environment or default to 5000
    app.run(host='0.0.0.0', port=port, debug=True)
    

# # Integrate user input
# if __name__ == "__main__":
#     user_input = input("Enter a DNA sequence: ")

#     try:
#         result = translate_dna(user_input)
#         print(f'Translated Protein Sequence: {result}')
#     except ValueError as e:
#         print(e)                              