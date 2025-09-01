# Magnolia

[![PyPI version](https://badge.fury.io/py/magnolia.svg)](https://badge.fury.io/py/magnolia)
[![Python](https://img.shields.io/pypi/pyversions/magnolia.svg)](https://pypi.org/project/magnolia/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Downloads](https://pepy.tech/badge/magnolia)](https://pepy.tech/project/magnolia)

Magnolia is a powerful Python package designed for comprehensive post-processing and analysis of LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator) molecular dynamics simulation data. It provides a streamlined interface for handling simulation outputs with robust data processing capabilities.

## Features

### Core Capabilities
- 📊 Parse and analyze LAMMPS log files with high performance
- 📈 Extract thermodynamic data into pandas DataFrames
- 🔄 Support for serial analysis with flexible data selection
- 🎯 Zero-referencing capabilities for energy and other parameters
- ⚡ Built-in timestep detection

### Advanced Features
- 📌 Automatic data validation and error checking
- 🔍 Comprehensive data analysis tools
- 📉 Built-in plotting and visualization functions
- 💾 Export capabilities to various formats
- 🔧 Extensible architecture for custom analysis

## Installation

You can install Magnolia using pip:

```bash
pip install magnolia
```

For development installation, clone the repository and install in editable mode:

```bash
git clone https://github.com/yourusername/magnolia.git
cd magnolia
pip install -e .
```

## Usage

Here's a simple example of how to use Magnolia to parse a LAMMPS log file:

```python
from magnolia import thermo_panda

# Parse a specific serial from the log file
data = thermo_panda(
    logfile="path/to/your/logfile",
    serial=1,  # Can be int, list, or "start:end" format
    zero_ref="energy"  # Optional zero-referencing
)

# Access your data as a pandas DataFrame
print(data.head())

# Basic analysis
average_energy = data['TotEng'].mean()
temperature_profile = data['Temp'].plot(title='Temperature Evolution')

# Export processed data
data.to_csv('processed_results.csv')
```

For more complex examples and detailed usage patterns, visit our [Examples Gallery](link_to_examples).

## Documentation

### Quick Start
For a quick introduction to Magnolia's features, check out our [Quick Start Guide](link_to_quickstart).

### API Reference
Complete API documentation is available at our [Documentation Portal](link_to_documentation).

### Examples
Find practical examples and tutorials in our [Examples Gallery](link_to_examples).

## Requirements

### Minimum Requirements
- Python 3.6+
- pandas >= 1.0.0
- numpy >= 1.18.0

### Optional Dependencies
- matplotlib >= 3.3.0 (for plotting capabilities)
- scipy >= 1.5.0 (for advanced analysis)
- pytest >= 6.0.0 (for running tests)

## Contributing

We welcome contributions to Magnolia! Here's how you can help:

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

Please ensure your PR adheres to our:
- [Code of Conduct](CODE_OF_CONDUCT.md)
- [Contributing Guidelines](CONTRIBUTING.md)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

- **Shihab Ahmed** - *Initial work*

## Citation

If you use Magnolia in your research, please cite:

```bibtex
@software{magnolia2024ahmed,
  author = {Ahmed, Shihab},
  title = {Magnolia: A post processing module for LAMMPS-MD simulation},
  year = {2024},
  url = {https://github.com/yourusername/magnolia}
}
```

## Acknowledgments

- LAMMPS development team
- Contributors and users of Magnolia
