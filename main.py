
#%%
import sys

sys.path.append("/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue")
from utils import *


class MutationsDataProcessor:
    def __init__(self, drug, genes_file, mutations_file):
        self.drug = drug
        self.genes_file = genes_file
        self.mutations_file = mutations_file
        self.genes = None

    def process_input_data(self):
        self.genes = get_genes_of_interest_from_genes_file(self.genes_file, self.drug)

        mutations_separator = get_separator(self.mutations_file)
        mutations = pd.read_csv(self.mutations_file, sep=mutations_separator, header=0).reset_index()

        check_mutations_file_columns(mutations)

        gene_column = next((col for col in mutations.columns if 'gene' in col.lower()), None)
        mutation_column = next((col for col in mutations.columns if 'mutation' in col.lower()), None)

        if gene_column:
            filtered_mutations = mutations[mutations[gene_column].isin(self.genes)].copy()
            filtered_mutations.loc[:, 'GENE_MUT'] = [
                f"{row[gene_column]}@{row[mutation_column]}"
                for _, row in filtered_mutations.iterrows()
            ]

            filtered_mutations['IS_SYNONYMOUS'] = [
                row[mutation_column][0] == row[mutation_column][-1]
                for _, row in filtered_mutations.iterrows()
            ]

            # Remove synonymous entries
            filtered_mutations = filtered_mutations[~filtered_mutations.IS_SYNONYMOUS]

            return filtered_mutations
        else:
            raise ValueError("Missing 'gene' column in the mutations file")


class PhenotypesDataProcessor:
    """
    A class to process phenotypes data.

    Attributes:
        drug (str): The drug of interest.
        genes_file (str): The file containing genes information.
        phenotypes_file (str): The file containing phenotypes data.
        phenotype_quality (str): The quality of phenotypes to filter ('HIGH', 'MEDIUM', or 'LOW').
        genes (list): List of genes of interest.
    """

    def __init__(self, drug, genes_file, phenotypes_file, phenotype_quality="MEDIUM"):
        """
        Initializes PhenotypesDataProcessor with the provided parameters.

        Args:
            drug (str): The drug of interest.
            genes_file (str): The file containing genes information.
            phenotypes_file (str): The file containing phenotypes data.
            phenotype_quality (str, optional): The quality of phenotypes to filter. Defaults to "MEDIUM".
                Possible values are "HIGH", "MEDIUM", or "LOW".
        """
        self.drug = drug
        self.genes_file = genes_file
        self.phenotypes_file = phenotypes_file
        self.phenotype_quality = phenotype_quality.upper()  # Convert to uppercase to handle case-insensitivity
        self.genes = None

    def process_input_data(self):
        """
        Processes the input data to filter phenotypes based on drug and phenotype quality.

        Returns:
            pd.DataFrame: DataFrame containing filtered phenotypes data.
        """
        self.genes = get_genes_of_interest_from_genes_file(self.genes_file, self.drug)

        phenotypes_separator = get_separator(self.phenotypes_file)
        phenotypes = pd.read_csv(self.phenotypes_file, sep=phenotypes_separator).reset_index()
        # filter for relevant drug
        # Check if any column contains the string "drug" (case-insensitive)
        drug_column = next((col for col in phenotypes.columns if 'drug' in col.lower()), None)

        # If a matching column is found, use it for filtering
        if drug_column:
            filtered_phenotypes = phenotypes[phenotypes[drug_column] == self.drug][["ena_run", "drug", "phenotype",
                                                                                    "phenotype_quality"]].copy()
        else:
            # Handle case when no matching column is found
            print("No 'drug' or 'DRUG' column found in the phenotypes file.")
            sys.exit()

        # Filter for 'R' and 'S' phenotypes only (excludes 'U' and NA)
        # Check if any column contains the string "phenotype" (case-insensitive)
        phenotype_column = next((col for col in filtered_phenotypes.columns if 'phenotype' in col.lower()), None)
        # If a matching column is found, use it for filtering
        if phenotype_column:
            filtered_phenotypes = filtered_phenotypes[filtered_phenotypes[phenotype_column].isin(["R", "S"])].copy()

            # Apply additional filters based on phenotype_quality
            if self.phenotype_quality == "HIGH":
                filtered_phenotypes = filtered_phenotypes[filtered_phenotypes["phenotype_quality"] == "HIGH"]
            elif self.phenotype_quality == "MEDIUM":
                filtered_phenotypes = filtered_phenotypes[
                    filtered_phenotypes["phenotype_quality"].isin(["HIGH", "MEDIUM"])]
            elif self.phenotype_quality == "LOW":
                filtered_phenotypes = filtered_phenotypes[
                    filtered_phenotypes["phenotype_quality"].isin(["HIGH", "MEDIUM", "LOW"])]
            else:
                print("Invalid phenotype_quality specified.")
                sys.exit()

            filtered_phenotypes = filtered_phenotypes.groupby("ena_run").apply(filter_multiple_phenos).reset_index(
                drop=True)
            return filtered_phenotypes
        else:
            print("No 'drug' or 'DRUG' column found in the phenotypes file.")
            sys.exit()


class CreateInputDataSummaryTable:
    def __init__(self, mutations_df, filtered_phenotypes_df, frs_threshold, drug):
        if not 0 <= frs_threshold <= 1:
            raise ValueError("frs_threshold must be a number between 0 and 1")
        self.mutations_df = mutations_df
        self.phenotypes_df = filtered_phenotypes_df
        self.frs_threshold = frs_threshold
        self.drug = drug

    def create_summary_table(self):
        self.all_data = pd.merge(self.mutations_df, self.phenotypes_df, on='ena_run', how='inner')
        self.all_data['gene'].fillna('None', inplace=True)

        df = RSIsolateTable(self.all_data, self.all_data.gene.unique())
        df1 = RSIsolateTable(self.all_data[self.all_data.frs < self.frs_threshold], self.all_data.gene.unique())
        df2 = RSVariantTable(self.all_data, self.all_data.gene.unique())
        df3 = RSVariantTable(self.all_data[self.all_data.frs < self.frs_threshold], self.all_data.gene.unique())
        df = pd.concat([df, df1, df2, df3], axis=1)

        # Modify column names based on frs_threshold
        threshold_str = f" (frs < {self.frs_threshold})"
        df.columns = pd.MultiIndex.from_tuples(
            zip(
                [
                    "All isolates",
                    "",
                    "",
                    f"Minor alleles isolates{threshold_str}",
                    "",
                    "",
                    "All variants",
                    "",
                    "",
                    f"Minor alleles variants{threshold_str}",
                    "",
                    "",
                ],
                df.columns,
            )
        )

        # Save the table to the working directory with the desired filename format
        filename = f"{self.drug}_frs{self.frs_threshold}_summary_counts.csv"
        df.to_csv(filename, index=False)
        print(f"Summary table saved to {filename}")

        return df


#     ### TODO: cat builder note where frs filtering and minor stuff can be added, test code.
#     # TODO: fix argparse code when fin
#     # TODO: tidy up utils.py and remove stuff not needed that DA had there

class BuildResCatalogueProcessor:
    def __init__(self, mutations_df, filtered_phenotypes_df, genbank_file, drug, frs_thresholds, catalogue_version,
                 piezo_wildcards=None, out_dir=None, hardcoded=None, merge_type='both'):
        """
        BuildResCatalogue class for constructing a resistance catalogue.

        Parameters:
            mutations_df (DataFrame): DataFrame containing mutations data.
            filtered_phenotypes_df (DataFrame): DataFrame containing filtered phenotypes data.
            genbank_file (str): Path to the GenBank file.
            drug (str): Name of the drug.
            frs_thresholds (list): List of frs thresholds.
            catalogue_version (str): Version of the catalogue.
            piezo_wildcards (dict, optional): Dictionary containing piezo wildcards. Defaults to None.
            out_dir (str, optional): Output directory path. Defaults to None.
            hardcoded (dict, optional): Dictionary containing hardcoded variant classifications. Defaults to None.
            merge_type (str, optional): Choose which samples from mutations_df and filtered_phenotypes_df to include.
                i.e. if you want to filter one df by another to choose a subset of samples.
                Choose 'mutations' (for ensuring mutations df ena_runs (samples) are in resulting DataFrame),
                'phenotypes' (for ensuring phenotypes df ena_runs (samples) are in resulting DataFrame),
                or 'both' (for merging both mutations and phenotypes and ensuring all ena_runs (samples)
                present in resulting DataFrame).
                Defaults to 'both'.
        """
        self.catalogue = None
        self.mutations_df = mutations_df
        self.filtered_phenotypes_df = filtered_phenotypes_df
        self.genbank_file = genbank_file
        self.drug = drug
        self.frs_thresholds = frs_thresholds
        self.catalogue_version = catalogue_version
        self.piezo_wildcards = piezo_wildcards
        if out_dir is None:
            self.out_dir = os.getcwd()  # Set to current working directory
        else:
            self.out_dir = out_dir
        self.hardcoded = hardcoded
        self.merge_type = merge_type

        # Ensure frs_thresholds are within the valid range
        for threshold in self.frs_thresholds:
            if not 0 <= threshold <= 1:
                raise ValueError("frs_threshold must be a number between 0 and 1")

        # Merge mutations and phenotypes dfs based on the merge_type parameter
        if self.merge_type == 'mutations':
            # Ensure all ena_runs from mutations_df are present in the resulting DataFrame
            self.all_data = pd.merge(self.mutations_df, self.filtered_phenotypes_df['ena_run'], on=["ena_run"], how="left")
            # Check if all ena_runs are present in the resulting all_data DataFrame
            if len(self.all_data) != len(self.mutations_df):
                raise ValueError("Not all ena_runs from mutations_df are present in the resulting all_data DataFrame."
                                 "Check that all ena_runs are represented in phenotypes_df")
        elif self.merge_type == 'phenotypes':
            # Ensure all ena_runs from filtered_phenotypes_df are present in the resulting DataFrame
            self.all_data = pd.merge(self.mutations_df['ena_run'], self.filtered_phenotypes_df, on=["ena_run"], how="right")
            # Check if all ena_runs are present in the resulting all_data DataFrame
            if len(self.all_data) != len(self.filtered_phenotypes_df):
                raise ValueError("Not all ena_runs from phenotypes_df are present in the resulting all_data DataFrame."
                                 "Check that all ena_runs are represented in mutations_df")
        elif self.merge_type == 'both':
            # Merge both mutations_df and filtered_phenotypes_df on "ena_run" column
            self.all_data = pd.merge(self.mutations_df, self.filtered_phenotypes_df, on=["ena_run"], how="inner")
        else:
            raise ValueError("Invalid merge_type. Choose from 'mutations', 'phenotypes', or 'both'.")
        unique_genes = self.all_data[~self.all_data.GENE_MUT.isna()].GENE_MUT.unique()
        #print(f"list of unique mutations: {unique_genes}")
        print(f"Length of unique mutations list: {len(unique_genes)}")
        print(f"All ready to make catalogue(s)...")

        if self.piezo_wildcards:
            if drug == "isoniazid":
                print(f"Adding Piezo wildcards for {drug}")
                self.piezo_wildcards = {
                    "inhA@*=": {"pred": "S", "evid": {}},
                    "inhA@-*_indel": {"pred": "U", "evid": {}},
                    "inhA@*_indel": {"pred": "U", "evid": {}},
                    "inhA@-*?": {"pred": "U", "evid": {}},
                    "inhA@*?": {"pred": "U", "evid": {}},
                    "inhA@del_0.0": {"pred": "U", "evid": {}},
                    "katG@*=": {"pred": "S", "evid": {}},
                    "katG@-*_indel": {"pred": "U", "evid": {}},
                    "katG@*_indel": {"pred": "U", "evid": {}},
                    "katG@-*?": {"pred": "U", "evid": {}},
                    "katG@*?": {"pred": "U", "evid": {}},
                    "katG@del_0.0": {"pred": "U", "evid": {}},
                }


    def build_catalogue(self):
        # Create output directory if it doesn't exist
        if self.out_dir and not os.path.exists(self.out_dir):
            print("Creating output directory...")
            os.makedirs(self.out_dir)

        for threshold in self.frs_thresholds:
            # Create self.catalogue_name from the input drug and frs_threshold string
            self.catalogue_name = f"{self.drug}_FRS{threshold}_catalogue"
            print(f"Building {self.catalogue_name}")
            output_csv_file = os.path.join(self.out_dir, f"{self.catalogue_name}.csv")

            # apply frs threshold to mutations
            self.all_data_frs_filtered = self.all_data[(self.all_data.frs >= threshold)].copy()
            n_variants = self.all_data_frs_filtered.shape[0]
            print(f"the number of variants for FRS threshold {threshold} is {n_variants}")

            self.S, self.R, self.U = [], [], []

            # hardcode variant classifications - often helps to seed with phylogenetic mutations
            if self.hardcoded:
                for k, v in self.hardcoded.items():
                    if v == 'S':
                        self.S.append({'mut': k, 'evid': {}})
            #print(f"length of self.S: {len(self.S)}")

            self.run = True
            # Run build_S_arr as long as self.run is True
            while self.run:
                build_S_arr(self)

            # Run mop_up if needed
            #if not self.run:
            print(f"All SOLOs found for S")
            mop_up(self)

            print(f"All done for finding S and R samples...")

            self.catalogue = construct_catalogue(self)
            if len(self.U) == 0:
                prediction_value_string = "RS"
            elif len(self.U) > 0:
                prediction_value_string = "RUS"
            # Output catalogue in piezo format
            self.piezo_catalogue = return_piezo(self, self.genbank_file, self.catalogue_name,
                                                               self.catalogue_version, self.drug,
                                                               self.piezo_wildcards, values=prediction_value_string).to_csv(output_csv_file)

    def return_all_data(self):
        return self.all_data



class PredictResistance:
    def __init__(self, mutations_df, phenotypes_df, catalogue_file, catalogue_name,
                 drug, frs_threshold, U_to_R=False, U_to_S=False,
                 merge_type="both", out_dir=None, Print=True):

        """
        Predict resistance given a df and resistance catalogue (piezo format)
        """
        self.mutations_df = mutations_df
        self.phenotypes_df = phenotypes_df
        self.catalogue_file = catalogue_file
        self.catalogue_name = catalogue_name
        self.drug = drug
        self.frs_threshold = frs_threshold
        self.U_to_R = U_to_R
        self.U_to_S = U_to_S
        self.merge_type = merge_type
        if out_dir is None:
            self.out_dir = os.getcwd()  # Set to current working directory
        else:
            self.out_dir = out_dir
        self.Print = Print


        self.out_csv_matrix = f"{self.catalogue_name}_{self.drug}_FRS{self.frs_threshold:.2f}_predictions_matrix.csv"
        self.out_csv_matrix_file = os.path.join(self.out_dir, f"{self.out_csv_matrix}")

        self.out_csv_predictions = f"{self.catalogue_name}_{self.drug}_FRS{self.frs_threshold:.2f}_predictions.csv"
        self.out_csv_predictions_file = os.path.join(self.out_dir, f"{self.out_csv_predictions}")

        self.out_text_error = f"{self.catalogue_name}_{self.drug}_FRS{self.frs_threshold:.2f}_predictions_errors.txt"
        self.out_text_error_file = os.path.join(self.out_dir, f"{self.out_text_error}")

        self.out_text_summary = f"{self.catalogue_name}_{self.drug}_FRS{self.frs_threshold:.2f}_predictions_performance_summary.txt"
        self.out_text_summary_file = os.path.join(self.out_dir, f"{self.out_text_summary}")

        # Merge mutations and phenotypes dfs based on the merge_type parameter
        if self.merge_type == 'mutations':
            # Ensure all ena_runs from mutations_df are present in the resulting DataFrame
            self.all_data = pd.merge(self.mutations_df, self.phenotypes_df['ena_run'], on=["ena_run"], how="left")
            # Check if all ena_runs are present in the resulting all_data DataFrame
            if len(self.all_data) != len(self.mutations_df):
                raise ValueError("Not all ena_runs from mutations_df are present in the resulting all_data DataFrame."
                                 "Check that all ena_runs are represented in phenotypes_df")
        elif self.merge_type == 'phenotypes':
            # Ensure all ena_runs from filtered_phenotypes_df are present in the resulting DataFrame
            self.all_data = pd.merge(self.mutations_df['ena_run'], self.phenotypes_df, on=["ena_run"], how="right")
            # Check if all ena_runs are present in the resulting all_data DataFrame
            if len(self.all_data) != len(self.phenotypes_df):
                raise ValueError("Not all ena_runs from phenotypes_df are present in the resulting all_data DataFrame."
                                 "Check that all ena_runs are represented in mutations_df")
        elif self.merge_type == 'both' or None:
            # Merge both mutations_df and filtered_phenotypes_df on "ena_run" column
            self.all_data = pd.merge(self.mutations_df, self.phenotypes_df, on=["ena_run"], how="inner")
        else:
            raise ValueError("Invalid merge_type. Choose from 'mutations', 'phenotypes', or 'both'.")



        # Filter data by frs if provided
        if self.frs_threshold:
            if not 0 <= self.frs_threshold <= 1:
                raise ValueError("frs_threshold must be a number between 0 and 1")
            self.all_data = self.all_data[(self.all_data.frs >= frs_threshold)].copy()
        else:
            pass


    def predict_resistance(self):
        cm = PiezoPredict(self.all_data, self.catalogue_file, self.out_csv_matrix_file,
                          self.out_text_error_file, self.out_text_summary_file, self.out_csv_predictions_file, self.drug)[0]
        return cm



# sens = 100 * (cm[0][0] / (cm[0][0] + cm[0][1]))
# spec = 100 * (cm[1][1] / (cm[1][1] + cm[1][0]))
# cov = 100 * (len(samples) - cm[0][2] - cm[1][2]) / len(samples)
# df_cm = pd.DataFrame(cm, index=["R", "S"], columns=["R", "S", "U"])
#
# plt.figure(figsize=(8, 5))
# sns.heatmap(
#     df_cm, annot=True, cbar=False, fmt="g", cmap="Greens", annot_kws={"fontsize": 24}
# )
# plt.gca().invert_yaxis()
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.xlabel("Predicted", fontsize=14)
# plt.ylabel("True", fontsize=14)
# plt.show()


#
#
# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Process some integers.")
#     parser.add_argument("drug", type=str, help="Enter the drug name of interest (in format relevant to "
#                                                "genes_file 'drug' column and 'drug' column in mutations_file)")
#     parser.add_argument("genes_file", type=str, help="path to the file containing a list of genes "
#                                                      "to be investigated. File columns required: [drug] [gene]")
#     parser.add_argument("mutations_file", type=str, help="path to the mutations file")
#     #parser.add_argument("genomes_file", type=str, help="path to the genomes file")
#     parser.add_argument("phenotypes_file", type=str, help="Path to the phenotypes file. Phenotypes file. "
#                                                           "Must have the following columns: ena_run, DRUG, "
#                                                           "phenotype, phenotype_QUALITY")
#     parser.add_argument("phenotype_quality", type=check_phenotype_quality_type, help="Phenotype quality: "
#                                                                                "HIGH, MEDIUM, or LOW")
#     parser.add_argument("genbank_ref", type=str, help="genbank reference")
#     parser.add_argument("catalogue_name", type=str, help="name of the catalogue")
#     parser.add_argument("catalogue_version", type=str, help="version of the catalogue")
#     parser.add_argument("piezo_wildcards", type=str, help="Piezo wildcards")
#     parser.add_argument("output_csv_file", type=str, help="path to the output CSV file")
#
#     args = parser.parse_args()
#     processor = DataProcessor(args.mutations_file, args.genomes_file, args.phenotypes_file, args.genbank_ref, args.catalogue_name, args.catalogue_version, args.drug, args.piezo_wildcards, args.genes_file, args.output_csv_file)
#     processor.process_input_data()

#%%%
mutations_processor = MutationsDataProcessor("isoniazid",
                                             "/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/gene_panel_20240125.tsv",
                                             "/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/MUTATIONS_training_katG_inhA_ahpC_subset.csv")

mutations = mutations_processor.process_input_data()

phenotypes_processor = PhenotypesDataProcessor("isoniazid",
                                               "/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/gene_panel_20240125.tsv",
                                               "/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/training_data_phenotypes_20240125.tsv",
                                               "MEDIUM")

phenotypes = phenotypes_processor.process_input_data()

summarytable_processor = CreateInputDataSummaryTable(mutations, phenotypes, 1, "isoniazid")

summarytable = summarytable_processor.create_summary_table()
#%%%
build_cat_processor = BuildResCatalogueProcessor(mutations, phenotypes, "NC00962.3", "isoniazid",
                                                 [1.0, 0.8, 0.6], "1.1",
                                                 out_dir= "/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results",
                                                 hardcoded={'katG@R463L':'S'},
                                                 piezo_wildcards=True)

all_data = build_cat_processor.return_all_data()
build_cat_processor.build_catalogue()
#%%%
#training FRS 1.0 cat

for prediction_frs_threshold in [1.0, 0.8, 0.6]:
    prediction_processor = PredictResistance(mutations, phenotypes,
                                             "/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results/isoniazid_FRS1.0_catalogue.csv",
                                             catalogue_name="training_FRS1.0_catalogue",
                                             drug="isoniazid",
                                             frs_threshold=prediction_frs_threshold,
                                             out_dir="/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results")

    predictions = prediction_processor.predict_resistance()
    #%%%
#training FRS 0.8 cat
for prediction_frs_threshold in [1.0, 0.8, 0.6]:
    prediction_processor = PredictResistance(mutations, phenotypes,
                                             "/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results/isoniazid_FRS0.8_catalogue.csv",
                                             catalogue_name="training_FRS0.8_catalogue",
                                             drug="isoniazid",
                                             frs_threshold=prediction_frs_threshold,
                                             out_dir="/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results")

    predictions = prediction_processor.predict_resistance()


#%%
#training FRS 0.6 cat
for prediction_frs_threshold in [1.0, 0.8, 0.6]:
    prediction_processor = PredictResistance(mutations, phenotypes,
                                             "/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results/isoniazid_FRS0.6_catalogue.csv",
                                             catalogue_name="training_FRS0.6_catalogue",
                                             drug="isoniazid",
                                             frs_threshold=prediction_frs_threshold,
                                             out_dir="/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results")

    predictions = prediction_processor.predict_resistance()



#%%%
#WHO cat
for prediction_frs_threshold in [1.0, 0.8, 0.6]:
    prediction_processor = PredictResistance(mutations, phenotypes,
                                             "/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/NC_000962.3_WHO-UCN-GTB-PCI-2021.7_v1.0_GARC1_RUS.csv",
                                             catalogue_name="WHO_catalogue",
                                             drug="INH",
                                             frs_threshold=prediction_frs_threshold,
                                             out_dir="/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results")

    predictions = prediction_processor.predict_resistance()
#%%%








# #%%
# print(build_cat_processor.all_data_frs_filtered[build_cat_processor.all_data_frs_filtered['GENE_MUT'] == "katG@Y155C"])
#
#
# no_S_mutations = build_cat_processor.all_data_frs_filtered[~build_cat_processor.all_data_frs_filtered.GENE_MUT.isin(i["mut"] for i in build_cat_processor.S)].copy()
# #%%
# print(len(no_S_mutations))
# solos = no_S_mutations.groupby("ena_run").filter(lambda x: len(x) == 1)
#
# #%%
# print(len(solos))
# print(solos[~no_S_mutations.GENE_MUT.isna()].GENE_MUT.unique())
#
# #%%%
# with open("/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results/tmp/S_no_mopfilt.txt", 'w') as summary_file:
#     # Write the summary content
#     summary_file.write(f"{build_cat_processor.S}")
#
# with open("/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results/tmp/R_no_mopfilt.txt", 'w') as summary_file:
# # Write the summary content
#     summary_file.write(f"{build_cat_processor.R}")
#
# with open("/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results/tmp/U_mopfilt.txt", 'w') as summary_file:
#     # Write the summary content
#     summary_file.write(f"{build_cat_processor.U}")
#
#
# build_cat_processor.all_data_frs_filtered.to_csv("/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results/tmp/all_data_frs_filtered_no_mopfilt.csv", index=False)
#
# no_S_mutations.to_csv("/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results/tmp/all_data_frs_filtered_no_S_mutations_no_mopfilt.csv", index=False)
#
# solos.to_csv("/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results/tmp/all_data_frs_filtered_no_S_mutations_filtered_one_sample_no_mopfilt.csv", index=False)
#
# #%%%
# print(build_cat_processor.all_data_frs_filtered[build_cat_processor.all_data_frs_filtered['GENE_MUT'] == "katG@Y155C"])
# #%%%
# import pandas as pd
# mop_solos = pd.read_csv("/Users/kmalone/macbook_m1_backup/github/cryptic-catalogue/data/results/tmp/all_data_frs_filtered_no_mopfilt.csv",
#                         header=0)
#
# def fisher_binary(solos, mut):
#     # Count occurrences of "R" and "S" phenotypes for the mutation and without the mutation
#     R_count = len(solos[(solos['phenotype'] == "R") & (solos['GENE_MUT'] == mut)])
#     S_count = len(solos[(solos['phenotype'] == "S") & (solos['GENE_MUT'] == mut)])
#     R_count_no_mut = len(solos[(solos['GENE_MUT'].isna()) & (solos['phenotype'] == "R")])
#     S_count_no_mut = len(solos[(solos['GENE_MUT'].isna()) & (solos['phenotype'] == "S")])
#
#     # Build contingency table: ((R count, S count), (background R count, background S count))
#     data = [[R_count, S_count], [R_count_no_mut, S_count_no_mut]]
#
#     # Calculate Fisher's exact test p-value
#     _, p_value = stats.fisher_exact(data)
#
#     # Determine the phenotype prediction based on the p-value and the uniqueness of the phenotype
#     if p_value < 0.05 or solos[solos['GENE_MUT'] == mut]['phenotype'].nunique() == 1:
#         # If variant frequency is 1 or p-value is significant, determine phenotype based on counts
#         pred = "R" if R_count > S_count else "S"
#     else:
#         pred = "U"  # If p-value is not significant or there is more than one phenotype, predict "U"
#
#     # Return the prediction along with evidence
#     return {
#         "pred": pred,
#         "evid": [[R_count, S_count], [R_count_no_mut, S_count_no_mut], [p_value, None]],
#     }



# mop_solos = build_cat_processor.all_data_frs_filtered[~build_cat_processor.all_data_frs_filtered.GENE_MUT.isin(i["mut"] for i in S)].copy()
#%%%
from scipy import stats
S, R, U = [], [], []

for mut in mop_solos[~mop_solos.GENE_MUT.isna()].GENE_MUT.unique():
        pheno = fisher_binary(mop_solos, mut)
        #print(f"prediction: {mut}\t{pheno}")
        if pheno["pred"] == "R":
            R.append({"mut": mut, "evid": pheno["evid"]})
        elif pheno["pred"] == "U":
            U.append({"mut": mut, "evid": pheno["evid"]})
#%%
