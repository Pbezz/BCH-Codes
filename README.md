# BCH Codes Project

This repository contains a Python implementation of general BCH Codes, developed as a mini project for the Information Theory course at the Mathematical Institute, University of Oxford. It is accompanied by a detailed report explaining the mathematical foundations, construction, encoding, and decoding of BCH Codes.

## Contents

- **BCH_Codes_report.pdf** – The full report detailing the theory and practical implementation of BCH Codes.
- **src/BCHCode.py** – A Python module implementing the `BCHCode` class for encoding and decoding messages.
- **examples/BCHCode_Tests.ipynb** – A Jupyter Notebook demonstrating tests and examples using the `BCHCode` class.

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
   ├── docs/
   │   └── BCH_Codes_report.pdf
   ├── tests/
   │   └── BCHCode_Tests.ipynb
   └── README.md
   ```

4. **Run the Examples:**
   ```bash
   jupyter notebook tests/BCHCode_Tests.ipynb
   ```

   Or import the module in your own Python script:
   ```python
   from src.BCHCode import BCHCode
   ```

## Usage

The `BCHCode`