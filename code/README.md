# Catalogue Generation and Analysis Tool

This tool generates a catalogue of genetic mutations associated with drug resistance in tuberculosis bacteria (Mycobacterium tuberculosis). It also performs analysis on the generated catalogue to predict drug resistance patterns and visualize the results.

## Requirements

- Python 3.x
- pandas
- matplotlib
- seaborn
- catalogue_builder (custom module)
- utils (custom module)

## Installation

1. Clone this repository to your local machine:  
   ```  
   git clone https://github.com/yourusername/catalogue-tool.git  
   cd catalogue-tool
   ```  
 2. Install the required Python packages using pip:
    ```
    pip install -r requirements.txt
    ```


## Usage

Run the main script with the required arguments to generate the catalogue and perform analysis:

```
python main.py \
mutations_file \
genomes_file \
phenotypes_file \
genbank_ref \
catalogue_name \
catalogue_version \
drug piezo_wildcards \
genes_file \
output_csv_file

```


Replace the arguments with the appropriate file paths and values.

## Arguments

- `mutations_file`: Path to the mutations file (CSV format).
- `genomes_file`: Path to the genomes file (CSV format).
- `phenotypes_file`: Path to the phenotypes file (CSV format).
- `genbank_ref`: Genbank reference.
- `catalogue_name`: Name of the catalogue.
- `catalogue_version`: Version of the catalogue.
- `drug`: Drug name.
- `piezo_wildcards`: Piezo wildcards.
- `genes_file`: Path to the file containing a list of genes.
- `output_csv_file`: Path to the output CSV file.

## Contributing

Contributions are welcome! If you find any issues or have suggestions for improvement, please open an issue or submit a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
