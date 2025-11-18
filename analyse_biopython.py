# analyse_biopython.py
"""
This script checks if Biopython is installed and prints its version. If not installed, it will print an error message.
"""

try:
    import Bio
    print(f"Biopython is installed. Version: {Bio.__version__}")
except ImportError:
    print("Biopython is not installed. Please install it using 'pip install biopython'.")
