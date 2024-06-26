import argparse
import os
import pandas as pd

# Create an argument parser
parser = argparse.ArgumentParser(description='Process *.variants.csv and *.mutations.csv files per sample to add FRS data')

parser.add_argument('input_filename', type=str, help='Input file containing a list of filenames')
parser.add_argument('variants_directory', type=str, help='Directory containing *.variants.csv files')
parser.add_argument('mutations_directory', type=str, help='Directory containing *.mutations.csv files')
parser.add_argument('output_directory', type=str, help='Output directory for the result files')

args = parser.parse_args()

input_filename = args.input_filename
variants_directory = args.variants_directory
mutations_directory = args.mutations_directory
output_directory = args.output_directory


# Get sample names from a list of filepaths to run per sample
def extract_sample_name(filename_or_path):
    # Use os.path.basename to extract the filename from a full path
    filename = os.path.basename(filename_or_path)
    # Split the filename and take the first part
    return filename.split('.variants.csv')[0]

# Process *.variants.csv file per sample
def process_variants_file(sample_name, variants_directory, mutations_directory):
    variants_filename = os.path.join(variants_directory, f"{sample_name}.variants.csv")
    df_variants = pd.read_csv(variants_filename, escapechar='\\')

    # Create a new data frame with selected columns and extract FRS metrics
    # Create a "merge" column
    # fillna(0) being used as there are intergenic variants present that have neither gene nor gene_position information
    df_new = df_variants[['uniqueid', 'gene', 'gene_position', 'vcf_evidence']].copy()
    df_new['merge'] = df_new['gene'].fillna(0).astype(str)  + '_' + df_new['gene_position'].fillna(0).astype(int).astype(str)
    df_new['FRS'] = df_new['vcf_evidence'].apply(lambda x: float(x.split('"FRS": ')[1].split(',')[0]) if '"FRS":' in x else None)

    # Read *.mutations.csv file per sample
    mutations_filename = os.path.join(mutations_directory, f"{sample_name}.mutations.csv")
    df_mutations = pd.read_csv(mutations_filename)

    # Create a "merge" column
    # fillna(0) being used as there are intergenic variants present that have neither gene nor gene_position information
    df_mutations['merge'] = (df_mutations['gene'].fillna(0).astype(str) + '_' +
                             df_mutations['gene_position'].fillna(0).astype(int).astype(str))

    # Merge data frames on the "merge" column
    df_result = pd.merge(df_new, df_mutations, on='merge', how='inner').drop_duplicates(subset='merge')
    columns_to_exclude = ['uniqueid_y', 'gene_y', 'gene_position_y']
    df_selected = df_result.drop(columns=columns_to_exclude)
    df_selected.rename(columns={'gene_x': 'gene', 'uniqueid_x': 'uniqueid', 'gene_position_x': 'gene_position'}, inplace=True)


    return df_selected

# Step 8: Output result to a file
def output_result(df_result, sample_name, output_directory):
    output_filename = os.path.join(output_directory, f"{sample_name}.mutations_frs.csv")
    df_result.to_csv(output_filename, index=False, header=True)

# Read the list of filenames from an input file
#input_filename = '/hps/nobackup/iqbal/kmalone/cryptic_catalogue/gnomonicus_output/variants_files.fofn'
#variants_directory = '/hps/nobackup/iqbal/kmalone/cryptic_catalogue/gnomonicus_output/variants'
#mutations_directory = '/hps/nobackup/iqbal/kmalone/cryptic_catalogue/gnomonicus_output/mutations'
#output_directory = '/hps/nobackup/iqbal/kmalone/cryptic_catalogue/gnomonicus_output/mutations_frs'


#with open(input_filename, 'r') as file:
#    file_list = [line.strip() for line in file]

# Loop to run on all samples
#for filename in file_list:
#    sample_name = extract_sample_name(filename)
#    result_df = process_variants_file(sample_name, variants_directory, mutations_directory)
#    output_result(result_df, sample_name, output_directory)

# Run per sample
sample_name = extract_sample_name(input_filename)
result_df = process_variants_file(sample_name, variants_directory, mutations_directory)
output_result(result_df, sample_name, output_directory)
