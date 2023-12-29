
# Mixed-integer programming for reduction of genome-scale metabolic networks

This repository contains a Python-based tool for metabolic network reduction using COBRApy, ZIMPL, and SCIP.

## Features

- Load metabolic network models in MAT or SBML format.
- Interactive user interface for selecting protected reactions and metabolites.
- Customization of parameters like delta, bigM, and degree of freedom (DoF).
- Functionality to add specific metabolic network constraints.
- Automated solution generation for minimum subnetworks considering user-defined constraints.

## Requirements

- Python 3.x
- COBRApy
- Tkinter
- NumPy
- Pandas
- SCIP Optimization Suite
- ZIMPL

## Installation

Clone the repository:

```bash
git clone https://github.com/haleminhduong/MetNetworkReducer.git
```

Install the required Python packages:

```bash
pip install cobra numpy pandas
```

Ensure SCIP and ZIMPL are installed and configured on your system.

## Usage

Run the main Python script from the program's root directory to launch the GUI:

```bash
python reduce.py
```

Follow the GUI instructions to load models, set parameters, and run the analysis.

## GUI Instructions

- **Load Model:** Choose a metabolic network model file (.mat or .xml).
- **Select Protected Reactions/Metabolites:** Choose reactions and metabolites to protect in the resulting network.
- **Set Parameters:** Input values for delta, bigM, and DoF.
- **Add Functionalities:** Define specific functionalities or constraints for the network.
- **Run Pipeline:** Execute the analysis to find minimum subnetworks.

## File Descriptions

- `reduce.py`: Main Python script with GUI and analysis logic.
- `zimpl_txt/`: Directory containing ZIMPL files for model constraints and parameters.

## Contributing

Contributions are welcome! Please read our contributing guidelines for details on how to submit pull requests, report issues, or request features.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgements

Thanks to all the contributors and users of the COBRApy, SCIP, and ZIMPL communities.
