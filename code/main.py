

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

from catalogue_builder import BuildCatalogue
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
        mutations = pd.read_csv(self.mutations_file, sep=mutations_separator).reset_index()

        check_mutations_file_columns(mutations)

        gene_column = next((col for col in mutations.columns if 'gene' in col.lower()), None)

        if gene_column:
            filtered_mutations = mutations[mutations[gene_column].isin(self.genes)]
            filtered_mutations['GENE_MUT'] = [
                f"{row['GENE']}@{row['MUTATION']}"
                for _, row in filtered_mutations.iterrows()
            ]

            filtered_mutations['IS_SYNONYMOUS'] = [
                row['MUTATION'][0] == row['MUTATION'][-1]
                for _, row in filtered_mutations.iterrows()
            ]

            # Remove synonymous entries
            filtered_mutations = filtered_mutations[~filtered_mutations.IS_SYNONYMOUS]

            return filtered_mutations
        else:
            raise ValueError("Missing 'gene' column in the mutations file")


class PhenotypesDataProcessor:
    def __init__(self, drug, genes_file, phenotypes_file, phenotype_quality):
        self.drug = drug
        self.genes_file = genes_file
        self.phenotypes_file = phenotypes_file
        self.phenotype_quality = phenotype_quality
        self.genes = None

    def process_input_data(self):
        self.genes = get_genes_of_interest_from_genes_file(self.genes_file, self.drug)

        phenotypes_separator = get_separator(self.phenotypes_file)
        phenotypes = pd.read_csv(self.phenotypes_file, sep=phenotypes_separator).reset_index()

        # filter for relevant drug
        # Check if any column contains the string "drug" (case-insensitive)
        drug_column = next((col for col in phenotypes.columns if 'drug' in col.lower() or 'drug' in col.upper()), None)

        # If a matching column is found, use it for filtering
        if drug_column:
            filtered_phenotypes = phenotypes[phenotypes[drug_column].isin(self.drug)][["UNIQUEID", "DRUG", "PHENOTYPE",
                                                                                       "PHENOTYPE_QUALITY"]]
        else:
            # Handle case when no matching column is found
            print("No 'drug' or 'DRUG' column found in the phenotypes file.")
            sys.exit()

        # Filter for 'R' and 'S' phenotypes only (excludes 'U' and NA)
        # Check if any column contains the string "phenotype" (case-insensitive)
        phenotype_column = next((col for col in filtered_phenotypes.columns if 'phenotype' in col.lower() or 'phenotype' in col.upper()), None)
        # If a matching column is found, use it for filtering
        if phenotype_column:
            filtered_phenotypes = filtered_phenotypes[filtered_phenotypes[phenotype_column].isin(["R", "S"])]
            filtered_phenotypes = filtered_phenotypes.groupby("UNIQUEID").apply(filter_multiple_phenos).reset_index(drop=True)
            return filtered_phenotypes

        else:
            # Handle case when no matching column is found
            print("No 'drug' or 'DRUG' column found in the phenotypes file.")
            sys.exit()




class CreateInputDataSummaryTable:
    def __init__(self, mutations_df, filtered_phenotypes_df, FRS_threshold):
        if not 0 <= FRS_threshold <= 1:
            raise ValueError("FRS_threshold must be a number between 0 and 1")
        self.mutations_df = mutations_df
        self.phenotypes_df = filtered_phenotypes_df
        self.FRS_threshold = FRS_threshold

    def create_summary_table(self):
        all_data = pd.merge(self.mutations_df, self.phenotypes_df, on='UNIQUEID', how='left')
        all_data['GENE'].fillna('None', inplace=True)

        df = RSIsolateTable(all_data, all_data.GENE.unique())
        df1 = RSIsolateTable(all_data[all_data.FRS < self.FRS_threshold], all_data.GENE.unique())
        df2 = RSVariantTable(all_data, all_data.GENE.unique())
        df3 = RSVariantTable(all_data[all_data.FRS < self.FRS_threshold], all_data.GENE.unique())
        df = pd.concat([df, df1, df2, df3], axis=1)

        df.columns = pd.MultiIndex.from_tuples(
            zip(
                [
                    "All",
                    "",
                    "",
                    "Minor alleles",
                    "",
                    "",
                    "All",
                    "",
                    "",
                    "Minor alleles",
                    "",
                    "",
                ],
                df.columns,
            )
        )
        return df





### TODO: cat builder class from stuff below, note where frs filtering and minor stuff can be added, test code.
        # TODO: fix argparse code when fin
        # TODO: tidy up utils.py and remove stuff not needed that DA had there

        #
        # samples = pd.DataFrame(samples)
        # samples['FRS'] = 1
        # cat_mutations = pd.DataFrame(cat_mutations)
        # cat_mutations['FRS'] = 1
        #
        # catalogue_90 = BuildCatalogue(samples, cat_mutations, 0.9)
        # data = catalogue_90.return_piezo(self.genbank_ref, self.catalogue_name, self.catalogue_version, self.drug, self.piezo_wildcards)
        # data.to_csv(self.output_csv_file)
        #
        # all_data = pd.merge(mutations, phenotypes, on=["UNIQUEID"], how="inner")
        # cm = PiezoPredict(all_data, self.output_csv_file, self.drug)["cm"]
        #
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





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument("drug", type=str, help="Enter the drug name of interest (in format relevant to "
                                               "genes_file 'drug' column and 'drug' column in mutations_file)")
    parser.add_argument("genes_file", type=str, help="path to the file containing a list of genes "
                                                     "to be investigated. File columns required: [drug] [gene]")
    parser.add_argument("mutations_file", type=str, help="path to the mutations file")
    #parser.add_argument("genomes_file", type=str, help="path to the genomes file")
    parser.add_argument("phenotypes_file", type=str, help="Path to the phenotypes file. Phenotypes file. "
                                                          "Must have the following columns: UNIQUEID, DRUG, "
                                                          "PHENOTYPE, PHENOTYPE_QUALITY")
    parser.add_argument("phenotype_quality", type=check_phenotype_quality_type, help="Phenotype quality: "
                                                                               "HIGH, MEDIUM, or LOW")
    parser.add_argument("genbank_ref", type=str, help="genbank reference")
    parser.add_argument("catalogue_name", type=str, help="name of the catalogue")
    parser.add_argument("catalogue_version", type=str, help="version of the catalogue")
    parser.add_argument("piezo_wildcards", type=str, help="Piezo wildcards")
    parser.add_argument("output_csv_file", type=str, help="path to the output CSV file")

    args = parser.parse_args()
    processor = DataProcessor(args.mutations_file, args.genomes_file, args.phenotypes_file, args.genbank_ref, args.catalogue_name, args.catalogue_version, args.drug, args.piezo_wildcards, args.genes_file, args.output_csv_file)
    processor.process_input_data()
