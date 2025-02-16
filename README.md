# BCH Codes Project

This repository contains a Python implementation of general BCH Codes, developed as a mini project for the Information Theory course at the Mathematical Institute, University of Oxford. It is accompanied by a detailed report explaining the mathematical foundations, construction, encoding, and decoding of BCH Codes.

## Contents

- **BCH_Codes_report.pdf** – The full report detailing the theory and practical implementation of BCH Codes.
- **BCHCode.py** – A Python module implementing the `BCHCode` class for encoding and decoding messages.
- **BCHCode_Tests.ipynb** – A Jupyter Notebook demonstrating tests and examples using the `BCHCode` class.
- **README.md** – This file.

## Getting Started

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/Pbezz/BCH-Codes
   ```
   
2. **Install Dependencies:**
   ```bash
   pip install sympy galois numpy jupyter
   ```

## Usage

The `BCHCode` class implements BCH codes over any finite field GF(q), where q is a prime number. Here's a basic example:

```python
from BCHCode import BCHCode

# Create a BCH code with parameters q=2, n=15, d=7, c=1
bch = BCHCode(2, 15, 7, 1)

# Encode a message
message = [1, 0, 1, 1, 0]
code = bch.encode(message)

# Decode (even with errors)
decoded = bch.decode(code)
```

For more examples and test cases, including text encoding/decoding, see `BCHCode_Tests.ipynb`.

## Project Overview
BCH Codes are cyclic error-correcting codes with strong error detection and correction capabilities. This project covers:

- The theoretical basis of BCH Codes, including finite fields  and polynomial codes.
- Detailed encoding and decoding algorithms.
- A complete Python implementation with interactive testing.

## Author

Developed by Pedro Costa as part of the Information Theory course at the University of Oxford.
