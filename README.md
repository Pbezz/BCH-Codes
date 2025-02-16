# BCH Codes Project

This repository contains a Python implementation of general BCH Codes, developed as a mini project for the Information Theory course at the Mathematical Institute, University of Oxford. It is accompanied by a detailed report explaining the mathematical foundations, construction, encoding, and decoding of BCH Codes.

## Contents

- **BCH_Codes_report.pdf** – The full report detailing the theory and practical implementation of BCH Codes.
- **src/BCHCode.py** – A Python module implementing the `BCHCode` class for encoding and decoding messages.
- **examples/BCHCode_Examples.ipynb** – A Jupyter Notebook demonstrating tests and examples using the `BCHCode` class.

## Getting Started

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/Pbezz/BCH-Codes
   cd BCH-Codes
   ```
   
2. **Install Dependencies:**
   ```bash
   pip install sympy galois numpy jupyter
   ```

3. **Project Structure:**
   Ensure your directory structure looks like this:
   ```
   BCH-Codes/
   ├── src/
   │   └── BCHCode.py
   ├── examples/
   │   └── BCHCode_Examples.ipynb
   ├── BCH_Codes_report.pdf
   └── README.md
   ```

4. **Run the Examples:**
   ```bash
   jupyter notebook examples/BCHCode_Examples.ipynb
   ```

   Or import the module in your own Python script:
   ```python
   from src.BCHCode import BCHCode
   ```

## Usage

The `BCHCode` class implements BCH codes over any finite field GF(q), where q is a prime number. Here's a basic example:

```python
# Create a BCH code with parameters q=2, n=15, d=7, c=1
bch = BCHCode(2, 15, 7, 1)

# Encode a message
message = [1, 0, 1, 1, 0]
code = bch.encode(message)

# Decode (even with errors)
decoded = bch.decode(code)
```

For more examples and test cases, including text encoding/decoding, see `BCHCode_Examples.ipynb`.

## Project Overview
BCH Codes are cyclic error-correcting codes with strong error detection and correction capabilities. This project covers:

- The theoretical basis of BCH Codes, including finite fields  and polynomial codes.
- Detailed encoding and decoding algorithms.
- A complete Python implementation with interactive testing.

## Author

Developed by Pedro Costa as part of the Information Theory course at the University of Oxford.
