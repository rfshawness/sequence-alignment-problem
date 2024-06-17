# Sequence Alignment Solutions

This project provides three distinct solutions to the sequence alignment problem:

1. Classic Dynamic Programming Solution
2. Memory Efficient Solution using Hirschberg's Algorithm
3. Miniature Version of the BLAST Algorithm (in progress)

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Directory Structure](#directory-structure)
- [Contributing](#contributing)
- [License](#license)

## Introduction

Sequence alignment is a fundamental problem in bioinformatics, essential for comparing DNA, RNA, or protein sequences. This repository offers implementations of three different approaches:

1. **Classic Dynamic Programming Solution:** A traditional approach using a matrix to compute the optimal alignment.
2. **Memory Efficient Solution using Hirschberg's Algorithm:** A space-efficient method that reduces the memory requirements of the dynamic programming solution.
3. **Miniature Version of the BLAST Algorithm:** An ongoing implementation of a simplified version of the BLAST algorithm, widely used for comparing sequence data.

## Installation

To get started, clone the repository and ensure you have Python 3 installed.

It's recommended to create a virtual environment to isolate your project and avoid conflicts with other packages. You can do this using the `venv` module that comes with Python.

```bash
python3 -m venv env
source env/bin/activate  # On Windows, use `env\Scripts\activate`
```

Once the virtual environment is activated, you can install the required packages. This project requires the `psutil` package.

```bash
pip3 install psutil
```

## Usage

Each solution is located in its respective directory. Within each directory, you will find a .py file for the Python implementation.  After navigating within the desired solution's directory, run the associated .py file with two arguments. The first argument is the relative path to your input.txt data file. The second argument being the relative path and name for the desired output files.

## Contributing

Contributions are welcome! If you have any improvements or new features to add, please follow these steps:

1. Fork the repository
2. Create a new branch (`git checkout -b feature-branch`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push to the branch (`git push origin feature-branch`)
5. Create a new Pull Request

## License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT). See the [LICENSE](LICENSE) file in this repository for details.