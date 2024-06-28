#%%
import re
import sys
import subprocess
import importlib.util
sys.path.append("/github/cryptic-catalogue")
from utils import *


class MutationsDataProcessor:
    def __init__(self, drug, genes_file, mutations_file, prepare_catomatic_input=False):
        self.drug = drug
        self.genes_file = genes_file
        self.mutations_file = mutations_file
        self.genes = None
        self.prepare_catomatic_input = prepare_catomatic_input

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

            if not self.prepare_catomatic_input:
                return filtered_mutations

            elif self.prepare_catomatic_input:
                catomatic_input = prepare_catomatic_mutations_input(filtered_mutations)
                return catomatic_input

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

    def __init__(self, drug, genes_file, phenotypes_file, phenotype_quality="MEDIUM", prepare_catomatic_input=False):
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
        self.prepare_catomatic_input =prepare_catomatic_input

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

            if not self.prepare_catomatic_input:
                return filtered_phenotypes

            elif self.prepare_catomatic_input:
                catomatic_input = prepare_catomatic_phenotypes_input(filtered_phenotypes)
                return catomatic_input
        else:
            print("No 'drug' or 'DRUG' column found in the phenotypes file.")
            sys.exit()


class CreateInputDataSummaryTable:
    def __init__(self, mutations_df, filtered_phenotypes_df, frs_threshold, drug, out_dir):
        if not 0 <= frs_threshold <= 1:
            raise ValueError("frs_threshold must be a number between 0 and 1")
        self.mutations_df = mutations_df
        self.phenotypes_df = filtered_phenotypes_df
        self.frs_threshold = frs_threshold
        self.drug = drug
        self.out_dir = out_dir

        print(f"Creating results directory for {drug}")
        self.results_dir = os.path.join(self.out_dir, f"{self.drug}")
        os.makedirs(self.results_dir, exist_ok=True)


    def create_summary_table(self):
        self.all_data = pd.merge(self.mutations_df, self.phenotypes_df, on='ena_run', how='inner')
        self.all_data['gene'].fillna('None', inplace=True)

        df = RSIsolateTable(self.all_data, self.all_data.gene.unique())
        df1 = RSIsolateTable(self.all_data[self.all_data.frs >= self.frs_threshold], self.all_data.gene.unique())
        df2 = RSVariantTable(self.all_data, self.all_data.gene.unique())
        df3 = RSVariantTable(self.all_data[self.all_data.frs >= self.frs_threshold], self.all_data.gene.unique())
        df = pd.concat([df, df1, df2, df3], axis=1)

        # Modify column names based on frs_threshold
        threshold_str = f" (frs >= {self.frs_threshold})"
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
        out_filepath = os.path.join(self.results_dir, filename)
        df.to_csv(out_filepath, index=False)
        print(f"Summary table saved to {filename}")

        return df


class BuildResCatalogueProcessor:
    def __init__(self, mutations_df, filtered_phenotypes_df, genbank_file, drug, frs_thresholds, catalogue_version,
                 run_OR =False, piezo_wildcards=None, out_dir=None, hardcoded_file=None, merge_type='both'):
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
        self.run_OR = run_OR
        self.piezo_wildcards = piezo_wildcards
        if out_dir is None:
            self.out_dir = os.getcwd()  # Set to current working directory
        else:
            self.out_dir = out_dir
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
            cwd = os.getcwd()
            wildcards_file = os.path.join(cwd, "data", "wildcards.json")
            with open(wildcards_file, 'r') as file:
                self.piezo_wildcards_data = json.load(file)
            print(f"Adding Piezo wildcards for {self.drug}")
            if self.drug in self.piezo_wildcards_data:
                self.piezo_wildcards = self.piezo_wildcards_data[self.drug]
                print(f"Setting Piezo wildcards for {self.drug}")
                self.prediction_value_string = "RUS"
            else:
                print(f"No Piezo wildcards found for {self.drug}")
                self.piezo_wildcards = None

        self.hardcoded_file = hardcoded_file
        if self.hardcoded_file is not None:
            # cwd = os.getcwd()
            # hardcoded_file = os.path.join(cwd, "data", "neutral_variants.json")
            with open(self.hardcoded_file, 'r') as file:
                self.hardcoded_data = json.load(file)
            if self.drug in self.hardcoded_data:
                self.hardcoded_data = self.hardcoded_data[self.drug]
                print(f"Found neutral variants for {self.drug}")
                self.prediction_value_string = "RUS"
            else:
                print(f"No neutral variants found for {self.drug}")
                self.hardcoded_file = None

        self.results_dir = os.path.join(self.out_dir, f"{self.drug}")
        os.makedirs(self.results_dir, exist_ok=True)


    def build_catalogue(self):
        # Create output directory if it doesn't exist
        if self.out_dir and not os.path.exists(self.out_dir):
            print("Creating output directory...")
            os.makedirs(self.out_dir)

        for threshold in self.frs_thresholds:
            # Create self.catalogue_name from the input drug and frs_threshold string
            self.catalogue_name = f"{self.drug}_FRS{threshold}_catalogue"
            print(f"Building {self.catalogue_name}")
            output_csv_file = os.path.join(self.results_dir, f"{self.catalogue_name}.csv")

            # apply frs threshold to mutations
            self.all_data_frs_filtered = self.all_data[(self.all_data.frs >= threshold)]
            n_variants = self.all_data_frs_filtered.shape[0]
            print(f"the number of variants for FRS threshold {threshold} is {n_variants}")

            self.S, self.R, self.U, self.to_remove = [], [], [], []

            # hardcode variant classifications - often helps to seed with phylogenetic mutations
            if self.hardcoded_data:
                for k, v in self.hardcoded_data.items():
                    if v == 'S':
                        self.S.append({'mut': k, 'evid': {}})

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
            if len(self.U) == 0 and not self.prediction_value_string:
                self. prediction_value_string = "RS"
            elif len(self.U) > 0:
                self.prediction_value_string = "RUS"

            print(f"Length of S: {len(self.S)}")
            print(f"Length of R: {len(self.R)}")
            print(f"Length of U: {len(self.U)}")
            # Output catalogue in piezo format
            self.piezo_catalogue = return_piezo(self, self.genbank_file, self.catalogue_name,
                                                               self.catalogue_version, self.drug,
                                                               self.piezo_wildcards, values=self.prediction_value_string).to_csv(output_csv_file)

    def return_all_data(self):
        return self.all_data



class BuildCatomaticCatalogueProcessor:
    def __init__(self, samples, mutations, genbank_ref, catalogue_name,
                 version, drug, wildcards=None, path_to_catomatic=None,
                 seeds=False, out_dir=False, test='Binomial', background=0.01, p=0.90,
                 FRS=1.00, save_to_piezo=False):
        self.samples = samples
        self.mutations = mutations
        self.path_to_catomatic =os.getcwd() if path_to_catomatic is None else path_to_catomatic
        self.out_dir = os.getcwd() if out_dir is None else out_dir
        self.test = test
        self.background = background
        print(f'background: {self.background}')
        self.p = p
        print(f"p: {self.p}")
        self.FRS = FRS
        self.genbank_ref = genbank_ref
        self.catalogue_name = catalogue_name
        self.version = version
        self.drug = drug
        self.seeds = seeds
        if self.seeds is not None:
            with open(self.seeds, 'r') as file:
                self.seeds_all = json.load(file)
        if self.drug in self.seeds_all:
            self.seeds = list(self.seeds_all[self.drug].keys())
            print(f"Found neutral variants for {self.drug}")
        else:
            print(f"No neutral variants found for {self.drug}")
        self.wildcards = wildcards
        if self.wildcards:
            with open(wildcards, 'r') as file:
                self.wildcards_data = json.load(file)
            print(f"Adding Piezo wildcards for {self.drug}")
            if self.drug in self.wildcards_data:
                self.wildcards = self.wildcards_data[self.drug]
                print(f"Setting Piezo wildcards for {self.drug}")
            else:
                print(f"No Piezo wildcards found for {self.drug}")
                self.wildcards = None
        self.save_to_piezo = save_to_piezo

        self.catalogue_builder = None

        # Check for the repository and import the module
        self.check_repo()
        self.import_catalogue_builder()

    def check_repo(self):
        if os.path.exists(self.path_to_catomatic):
            print(f"The directory {self.path_to_catomatic} exists. Proceeding with loading the module.")
        else:
            raise FileNotFoundError(f"The directory {self.path_to_catomatic} does not exist. Please clone the repository manually.")

    def import_catalogue_builder(self):
        module_path = os.path.join(self.path_to_catomatic, 'src', 'catomatic', 'CatalogueBuilder.py')
        if os.path.exists(module_path):
            print(f"The catomatic repo directory {self.path_to_catomatic} exists. Proceeding with loading the module.")
        else:
            raise FileNotFoundError(f"The file {module_path} does not exist. Please ensure the repository is correctly cloned.")

        spec = importlib.util.spec_from_file_location('CatalogueBuilder', module_path)
        catalogue_builder = importlib.util.module_from_spec(spec)
        sys.modules['CatalogueBuilder'] = catalogue_builder
        spec.loader.exec_module(catalogue_builder)
        self.catalogue_builder = catalogue_builder
        print("CatalogueBuilder module imported")

    def build_catalogue(self):
        if self.catalogue_builder is None:
            raise ImportError("CatalogueBuilder module is not imported.")

        # Get the BuildCatalogue function from catomatic
        BuildCatalogue = getattr(self.catalogue_builder, 'BuildCatalogue', None)
        if BuildCatalogue is None:
            raise AttributeError("The function 'BuildCatalogue' is not found in the CatalogueBuilder module.")

        # Call the BuildCatalogue function
        catalogue=BuildCatalogue(samples=self.samples, mutations=self.mutations, test=self.test, background=self.background,
                                 p=self.p, FRS=self.FRS)

        # Write out piezo catalogue
        if self.save_to_piezo:
            results_dir = os.path.join(self.out_dir, f"{self.drug}")
            os.makedirs(results_dir, exist_ok=True)
            outfile = os.path.join(results_dir, f"{self.drug}_catalogue_bg_{self.background}_p_{self.p}_FRS_{self.FRS}.csv")

            catalogue.to_piezo(genbank_ref=self.genbank_ref, catalogue_name=self.catalogue_name,
                               version=self.version, drug=drug, wildcards=self.wildcards, outfile=outfile)
            print("Catalogue has been converted to Piezo format and saved to:", outfile)

        return catalogue

    def create_final_catalogue(self):
        catalogue_files = []
        for root, _, files in os.walk(self.out_dir):
            for file in files:
                if file.endswith(f'_catalogue_bg_{self.background}_p_{self.p}_FRS_{self.FRS}.csv'):
                    catalogue_files.append(os.path.join(root, file))

        # Read the first file with header
        df = pd.read_csv(catalogue_files[0])
        df.reset_index(drop=True, inplace=True)

        # Read the rest of the files without header and concatenate
        for file in catalogue_files[1:]:
            temp_df = pd.read_csv(file, header=0)
            temp_df.reset_index(drop=True, inplace=True)
            df = pd.concat([df, temp_df], ignore_index=True)
            df.reset_index(drop=True, inplace=True)

        # Save the concatenated DataFrame to the output file
        outfile_catalogue = os.path.join(self.out_dir, f"catalogue_bg_{self.background}_p_{self.p}_FRS_{self.FRS}.csv")

        df.reset_index(drop=True, inplace=True)
        df.to_csv(outfile_catalogue, index=False)
        print(f"Concatenated file saved as: {outfile_catalogue}")



class PredictResistance:
    def __init__(self, mutations_df, phenotypes_df, catalogue_file, catalogue_name,
                 drug, frs_threshold, U_to_R=False, U_to_S=False,
                 merge_type="both", out_dir=None, outfile_prefix_string=None, Print=True):

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
        print(f"Output directory: {self.out_dir}")

        self.outfile_prefix_string = outfile_prefix_string
        print(f"Output file prefix: {self.outfile_prefix_string}")
        self.Print = Print

        self.results_dir = os.path.join(self.out_dir, f"{self.drug}")
        print(f"Results directory: {self.results_dir}")
        os.makedirs(self.results_dir, exist_ok=True)

        self.out_csv_matrix = f"{self.catalogue_name}_FRS_{self.frs_threshold}_{self.outfile_prefix_string}_predictions_matrix.csv"
        self.out_csv_matrix_file = os.path.join(self.results_dir, f"{self.out_csv_matrix}")

        self.out_csv_predictions = f"{self.catalogue_name}_FRS_{self.frs_threshold}_{self.outfile_prefix_string}_predictions.csv"
        self.out_csv_predictions_file = os.path.join(self.results_dir, f"{self.out_csv_predictions}")

        self.out_text_error = f"{self.catalogue_name}_FRS_{self.frs_threshold}_{self.outfile_prefix_string}_predictions_errors.txt"
        self.out_text_error_file = os.path.join(self.results_dir, f"{self.out_text_error}")

        if self.outfile_prefix_string is None:
            self.out_text_summary = f"FRS{self.frs_threshold:.2f}_predictions_performance_summary.txt"
        elif self.outfile_prefix_string is not None:
            self.out_text_summary = f"{self.catalogue_name}_FRS_{self.frs_threshold}_{self.outfile_prefix_string}_predictions_performance_summary.txt"
            self.out_text_summary_file = os.path.join(self.results_dir, f"{self.out_text_summary}")

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


        stat_summary = PiezoPredict(self.all_data, self.catalogue_file, self.out_csv_matrix_file,
                          self.out_text_error_file, self.out_text_summary_file, self.out_csv_predictions_file, self.drug,
                          self.U_to_R, self.U_to_S)

        background = self.outfile_prefix_string.split("_")[3]
        ci = self.outfile_prefix_string.split("_")[5]

        master_stats = f"{self.catalogue_name}_{self.drug}_prediction_stats.csv"
        master_stats_file = os.path.join(self.results_dir, f"{master_stats}")
        # Check if the master stats file exists
        if not os.path.exists(master_stats_file):
            # Open the master stats file in write mode and write the header
            with open(master_stats_file, 'w') as stats_file:
                stats_file.write("DRUG,BACKGROUND_RATE,CI,FRS_PREDICT,SENSITIVITY,SPECIFICITY,COVERAGE,#SS,#RR,#RS,#SR,#RU,#SU\n")

        # Open the output master stats file in append mode
        print(f"writing to master stats file: {master_stats_file}")
        with open(master_stats_file, 'a') as stats_file:
            # Write stats_summary to file
            stats_file.write(f"{self.drug},{background},{ci},{self.frs_threshold},{stat_summary}\n")




#%%
import traceback

drugs_mutations = {
    "amikacin":"/data/training/MUTATIONS_training_eis_rrs_subset.csv",
    "capreomycin":"/data/training/MUTATIONS_training_rrs_tlyA_subset.csv",
    "ciprofloxacin":"/data/training/MUTATIONS_training_gyrA_subset.csv",
    "delamanid":"/data/training/MUTATIONS_training_ddn_subset.csv",
    "ethambutol":"/data/training/MUTATIONS_training_embA_embB_subset.csv",
    "ethionamide":"/data/training/MUTATIONS_training_ethA_fabG1_inhA_subset.csv",
    "isoniazid":"/data/training/MUTATIONS_training_katG_inhA_ahpC_fabG1_subset.csv",
    "kanamycin":"/data/training/MUTATIONS_training_eis_rrs_subset.csv",
    "levofloxacin":"/data/training/MUTATIONS_training_gyrA_gyrB_subset.csv",
    "linezolid":"/data/training/MUTATIONS_training_rplC_subset.csv",
    "moxifloxacin":"/data/training/MUTATIONS_training_gyrA_gyrB_subset.csv",
    "ofloxacin":"/data/training/MUTATIONS_training_gyrA_subset.csv",
    "rifampicin":"/data/training/MUTATIONS_training_rpoB_subset.csv",
    "streptomycin":"/data/training/MUTATIONS_training_gid_rpsL_rrs_subset.csv"
}




with open('/data/drug_codes.json') as f:
    drug_codes = json.load(f)


# Iterate over each drug and its mutations file
for drug, mutations_file in drugs_mutations.items():
    try:
        mutations_processor = MutationsDataProcessor(drug=drug,
                                                     genes_file="/data/gene_panel_20240125.tsv",
                                                     mutations_file=mutations_file,
                                                     prepare_catomatic_input=True)

        mutations = mutations_processor.process_input_data()

        phenotypes_processor = PhenotypesDataProcessor(drug=drug,
                                                       genes_file="/data/gene_panel_20240125.tsv",
                                                       phenotypes_file="/data/training/training_data_phenotypes_20240125.tsv",
                                                       phenotype_quality="MEDIUM",
                                                       prepare_catomatic_input=True)

        phenotypes = phenotypes_processor.process_input_data()

        # Find the corresponding code for the drug
        if drug in drug_codes:
            drug_short = drug_codes[drug]
            print("Short code for", drug, "is", drug_short)
        else:
            print("Drug code not found")
            continue

        for this_background in [0.15, 0.175, 0.20, 0.225, 0.25]:
            for this_p in [0.90, 0.95, 0.99]:
                print(f"on background {this_background} and p {this_p}")
                build_catomatic_processor = BuildCatomaticCatalogueProcessor(
                    samples=phenotypes,
                    mutations=mutations,
                    out_dir='/data/training/results/catomatic_training',
                    path_to_catomatic="catomatic",
                    test='Binomial',
                    background=this_background,
                    p=this_p,
                    FRS=1,
                    genbank_ref="NC00962.3",
                    catalogue_name="catomatic_test",
                    version="1.1",
                    drug=drug,
                    seeds="/data/neutral_variants.json",
                    wildcards="/data/wildcards.json",
                    save_to_piezo=True)

                catomatic_catalogue = build_catomatic_processor.build_catalogue()


    except Exception as e:
        print(f"An error occurred while processing {drug}: {e}")
        traceback.print_exc()
        continue

# try:
#     print(f"Concatenating per-drug catalogues into one final catalogue")
#     build_catomatic_processor.create_final_catalogue()

# except Exception as e:
#     print(f"An error occurred while concatenating final catalogue: {e}")
#     traceback.print_exc()



#%%%
drugs_mutations = {
    "amikacin":"/data/training/MUTATIONS_training_eis_rrs_subset.csv",
    "capreomycin":"/data/training/MUTATIONS_training_rrs_tlyA_subset.csv",
    "ciprofloxacin":"/data/training/MUTATIONS_training_gyrA_subset.csv",
    "delamanid":"/data/training/MUTATIONS_training_ddn_subset.csv",
    "ethambutol":"/data/training/MUTATIONS_training_embA_embB_subset.csv",
    "ethionamide":"/data/training/MUTATIONS_training_ethA_fabG1_inhA_subset.csv",
    "isoniazid":"/data/training/MUTATIONS_training_katG_inhA_ahpC_fabG1_subset.csv",
    "kanamycin":"/data/training/MUTATIONS_training_eis_rrs_subset.csv",
    "levofloxacin":"/data/training/MUTATIONS_training_gyrA_gyrB_subset.csv",
    "linezolid":"/data/training/MUTATIONS_training_rplC_subset.csv",
    "moxifloxacin":"/data/training/MUTATIONS_training_gyrA_gyrB_subset.csv",
    "ofloxacin":"/data/training/MUTATIONS_training_gyrA_subset.csv",
    "rifampicin":"/data/training/MUTATIONS_training_rpoB_subset.csv",
    "streptomycin":"/data/training/MUTATIONS_training_gid_rpsL_rrs_subset.csv"
}



with open('/data/drug_codes.json') as f:
    drug_codes = json.load(f)

for drug, mutations_file in drugs_mutations.items():
    try:
        mutations_processor = MutationsDataProcessor(drug=drug,
                                                     genes_file="/data/gene_panel_20240125.tsv",
                                                     mutations_file=mutations_file,
                                                     prepare_catomatic_input=False)

        mutations = mutations_processor.process_input_data()

        phenotypes_processor = PhenotypesDataProcessor(drug=drug,
                                                       genes_file="/data/gene_panel_20240125.tsv",
                                                       phenotypes_file="/data/training/training_data_phenotypes_20240125.tsv",
                                                       phenotype_quality="MEDIUM",
                                                       prepare_catomatic_input=False)

        phenotypes = phenotypes_processor.process_input_data()

        # Find the corresponding code for the drug
        if drug in drug_codes:
            drug_short = drug_codes[drug]
            print("Short code for", drug, "is", drug_short)
        else:
            print("Drug code not found")
            continue

        results_directory = '/data/training/results/catomatic_training'
        this_drug_directory = os.path.join(results_directory, drug)
        these_files = []

        for root, _, files in os.walk(this_drug_directory):
            for file in files:
                if file.endswith('FRS_1.csv'):
                    these_files.append(os.path.join(root, file))

        for this_file in these_files:
            for prediction_frs_threshold in [1.0, 0.8, 0.6]:
                print(f"Running resistance prediction with threshold {prediction_frs_threshold} with catalogue {this_file}")
                prediction_processor = PredictResistance(mutations_df=mutations,
                                                         phenotypes_df=phenotypes,
                                                         catalogue_file=this_file,
                                                         catalogue_name="catalogue_v1.1_FRS1",
                                                         drug=drug,
                                                         frs_threshold=prediction_frs_threshold,
                                                         outfile_prefix_string=f"{os.path.splitext(os.path.basename(this_file))[0]}_binomial",
                                                         out_dir="/data/training/results/catomatic_training/")

                # predictions = prediction_processor.predict_resistance()
                # print(predictions)

        print(F"All done with {drug}")

    except Exception as e:
        print(f"An error occurred while processing {drug}: {e}")
        traceback.print_exc()
        continue

