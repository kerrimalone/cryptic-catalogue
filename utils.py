#%%
import csv

from piezo import catalogue
from sklearn.metrics import confusion_matrix
import piezo
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from scipy import stats
import json
import os
import pandas as pd
from pathlib import Path


def get_separator(filename):
    """
    Ascertain file separator
    """
    try:
        with open(filename, 'r') as file:
            dialect = csv.Sniffer().sniff(file.read(1024))
        return dialect.delimiter
    except csv.Error:
        raise ValueError(f"Unable to determine the separator for file '{filename}'")

def check_gene_file_columns(df):
    """
    Check that the columns in the dataframe from genes file are correct
    """
    required_columns = {"drug", "gene"}
    df.columns = df.columns.str.lower()  # Convert column names to lowercase
    if not set(df.columns) == required_columns:
        required_columns_lower = {col.lower() for col in required_columns}
        if not required_columns_lower <= set(df.columns):
            raise ValueError("Genes file must have the following columns: drug, gene")
    else:
        print("all good with genes file requirements")

def check_mutations_file_columns(df):
    """
   Check that the columns in the dataframe from mutations file are correct
   """
    required_columns = {"ena_run", "gene", "mutation"}
    df.columns = df.columns.str.lower()  # Convert column names to lowercase
    if not set(df.columns) == required_columns:
        required_columns_lower = {col.lower() for col in required_columns}
        if not required_columns_lower <= set(df.columns):
            raise ValueError("Hi genes file must have the following columns: drug, gene")
    else:
        print("all good with mutations file requirements")

def check_phenotypes_file_columns(df):
    """
   Check that the columns in the dataframe from phenotypes file are correct
   """
    required_columns = {"ena_run", "drug", "phenotype", "phenotype_quality"}
    df.columns = df.columns.str.lower()  # Convert column names to lowercase
    if not set(df.columns) == required_columns:
        raise ValueError("Phenotypes file must have the following columns: ena_run, drug, phenotype, phenotype_quality")
def get_genes_of_interest_from_genes_file(filename, drug):
    """
    Get list of genes by filtering genes list df by drug of interest
    """

    # if not isinstance(filename, str) or not os.path.exists(filename):
    #     raise ValueError("Filename must be a valid file path.")

    # Get separator or assume
    genes_separator = get_separator(filename)
    try:
        genes_df = pd.read_csv(filename, sep=genes_separator, header=0)
    except Exception as e:
        raise ValueError(f"Hi Error reading genes file: {e}")

    # Check columns
    check_gene_file_columns(genes_df)
    print("getting genes")

    genes_gene_column = next((col for col in genes_df.columns if 'gene' in col.lower()), None)
    genes_drug_column = next((col for col in genes_df.columns if 'drug' in col.lower()), None)

    if genes_gene_column and genes_drug_column:
        # Convert drug column to lowercase (or uppercase)
        genes_df[genes_drug_column] = genes_df[genes_drug_column].str.upper()

        # Convert user input to lowercase (or uppercase)
        drug = drug.upper()

        # Filter genes based on case-insensitive comparison
        filtered_genes = genes_df.loc[genes_df[genes_drug_column] == drug, genes_gene_column].tolist()
        return filtered_genes
    else:
        raise ValueError("Missing 'gene' or 'drug' column in the genes file")

def filter_multiple_phenos(group):
    """
    If a sample contains more than 1 phenotype,
    keep the resistant phenotype if there is one.
    (Developer: DA)
    """
    if len(group) > 1:
        if "R" in group["phenotype"].values:
            return group[group["phenotype"] == "R"].iloc[0:1]
        else:
            return group.iloc[0:1]
    else:
        return group

def RSIsolateTable(df, genes):
    """returns df of number of R and S isolates"""
    table = {}
    table["Total"] = {
        "R": df[df.phenotype == "R"].ena_run.nunique(),
        "S": df[df.phenotype == "S"].ena_run.nunique(),
        "Total": df.ena_run.nunique(),
    }
    for i in genes:
        d = df[df.gene == i]
        table[i] = {
            "R": d[d.phenotype == "R"].ena_run.nunique(),
            "S": d[d.phenotype == "S"].ena_run.nunique(),
            "Total": d[d.phenotype == "R"].ena_run.nunique()
                     + d[d.phenotype == "S"].ena_run.nunique(),
        }

    return pd.DataFrame.from_dict(table).T


def RSVariantTable(df, genes):
    """returns df of number of R and S variants"""
    table = {}
    table["Total"] = {
        "R": df[df.phenotype == "R"].ena_run.count(),
        "S": df[df.phenotype == "S"].ena_run.count(),
        "Total": df.ena_run.count(),
    }
    for i in genes:
        d = df[df.gene == i]
        table[i] = {
            "R": d[d.phenotype == "R"].ena_run.count(),
            "S": d[d.phenotype == "S"].ena_run.count(),
            "Total": d[d.phenotype == "R"].ena_run.count()
                     + d[d.phenotype == "S"].ena_run.count(),
        }

    return pd.DataFrame.from_dict(table).T




def CombinedDataTable(all):

    df = RSIsolateTable(all, all.GENE.unique())
    df1 = RSIsolateTable(all[all.FRS < 0.9], all.GENE.unique())
    df2 = RSVariantTable(all, all.GENE.unique())
    df3 = RSVariantTable(all[all.FRS < 0.9], all.GENE.unique())
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



def build_S_arr(self):
    # remove mutations predicted as susceptible from df (to potentially proffer additional, effective solos)
    mutations = self.all_data_frs_filtered[~self.all_data_frs_filtered.GENE_MUT.isin(i["mut"] for i in self.S)]

    # extract samples with only 1 mutation
    solos = mutations.groupby("ena_run").filter(lambda x: len(x) == 1)


    # method is jammed - end here.
    if len(solos) == 0:
        self.run = False

    mut_count = 0
    s_iters = 0
    S_count = 0
    U_count = 0
    R_count = 0
    #print(f"The number of unique mutations is {len(solos[~solos.GENE_MUT.isna()].GENE_MUT.unique())}")
    # for non WT or synonymous mutations
    for mut in solos[~solos.GENE_MUT.isna()].GENE_MUT.unique():
        mut_count +=1
        #print(f"this mutation: {mut}")
        # determine phenotype of mutation using Fisher's test
        # print(f"{s_iters}: {mut}")
        pheno = fisher_binary(solos, mut)
        #print(f"this prediction: {pheno}")
        if pheno["pred"] == "S":
            #print(f"adding to self.S")
            # if susceptible, add mutation to phenotype array
            self.S.append({"mut": mut, "evid": pheno["evid"]})
            s_iters += 1
        if pheno["pred"] == "S":
            S_count += 1
        if pheno["pred"] == "U":
            U_count += 1
        if pheno["pred"] == "R":
            R_count += 1


    if s_iters == 0:
        # if no susceptible solos (ie jammed) - move to mop up
        self.run = False

    # print(f" mutations processed: " + str(mut_count))
    # print(f"S count: " + str(S_count))
    # print(f"R count: " + str(R_count))
    # print(f"U count: " + str(U_count))

    # return self.run
def mop_up(self):
    # remove mutations predicted as susceptible from df (to potentially proffer additional, effective solos)
    no_S_mutations = self.all_data_frs_filtered[~self.all_data_frs_filtered.GENE_MUT.isin(i["mut"] for i in self.S)]

    # extract samples with only 1 mutation
    mop_solos = no_S_mutations.groupby("ena_run").filter(lambda x: len(x) == 1)
    #print(f"The number of unique R/U mutations is {len(mop_solos[~mop_solos.GENE_MUT.isna()].GENE_MUT.unique())}")

    # for non WT or synonymous mutations
    for mut in mop_solos[~mop_solos.GENE_MUT.isna()].GENE_MUT.unique():
        # determine phenotype of mutation using Fisher's test and add mutation to phenotype array (should be no S)
        pheno = fisher_binary(mop_solos, mut)
        if pheno["pred"] == "R":
            self.R.append({"mut": mut, "evid": pheno["evid"]})
        elif pheno["pred"] == "U":
            #print(F"length or U before adding {mut}: {len(self.U)}")
            self.U.append({"mut": mut, "evid": pheno["evid"]})
            #print(F"length or U after adding {mut}: {len(self.U)}")
            #print(f"{pheno['evid']}")


def fisher_binary(solos, mut):
    # Count occurrences of "R" and "S" phenotypes for the mutation and without the mutation
    R_count = len(solos[(solos['phenotype'] == "R") & (solos['GENE_MUT'] == mut)])
    S_count = len(solos[(solos['phenotype'] == "S") & (solos['GENE_MUT'] == mut)])
    R_count_no_mut = len(solos[(~solos['GENE_MUT'].isna()) & (solos['GENE_MUT'] != mut) & (solos['phenotype'] == "R")])
    S_count_no_mut = len(solos[(~solos['GENE_MUT'].isna()) & (solos['GENE_MUT'] != mut) & (solos['phenotype'] == "S")])

    # Build contingency table: ((R count, S count), (background R count, background S count))
    data = [[R_count, S_count], [R_count_no_mut, S_count_no_mut]]

    # Calculate Fisher's exact test p-value
    _, p_value = stats.fisher_exact(data)

    if p_value < 0.05 or solos[solos.GENE_MUT == mut].phenotype.nunique() == 1:
        # if variant frequency is 1 simply call the phenotype, otherwise call phenotype at 95% confidence
        if R_count > S_count:
            return {
                "pred": "R",
                "evid": [[R_count, S_count], [R_count_no_mut, S_count_no_mut], [p_value, _]],
            }
        else:
            return {
                "pred": "S",
                "evid": [[R_count, S_count], [R_count_no_mut, S_count_no_mut], [p_value, _]],
            }

    else:
        return {
            "pred": "U",
            "evid": [[R_count, S_count], [R_count_no_mut, S_count_no_mut], [p_value, _]],
        }


def construct_catalogue(self):
    catalogue = {}
    for i in self.S:
        catalogue[i["mut"]] = {"pred": "S", "evid": i["evid"]}
    for i in self.R:
        catalogue[i["mut"]] = {"pred": "R", "evid": i["evid"]}
    for i in self.U:
        catalogue[i["mut"]] = {"pred": "U", "evid": i["evid"]}

    return catalogue

def return_catalogue(self):
    return {
        mutation: {"phenotype": data["pred"]}
        for mutation, data in self.catalogue.items()
    }


def insert_wildcards(self, wildcards):
    if self.catalogue is not None and wildcards is not None:
        self.catalogue = {**self.catalogue, **wildcards}
    elif self.catalogue is None and wildcards is not None:
        self.catalogue = wildcards
    elif self.catalogue is not None and wildcards is None:
        # Do nothing or handle the case where wildcards is None
        pass
    else:
        # Both self.catalogue and wildcards are None, handle as needed
        pass

def return_piezo(
        self,
        genbank_file,
        catalogue_name,
        catalogue_version,
        drug,
        piezo_wildcards,
        grammar="GARC1",
        values="RUS",
):
    # insert piezo wildcards into catalogue object
    insert_wildcards(self,piezo_wildcards)

    piezo = (
        pd.DataFrame.from_dict(self.catalogue, orient="index")
        .reset_index()
        .rename(
            columns={"index": "MUTATION", "pred": "PREDICTION", "evid": "EVIDENCE", "p":"p_value"}
        )
    )
    piezo["GENBANK_REFERENCE"] = genbank_file
    piezo["CATALOGUE_NAME"] = catalogue_name
    piezo["CATALOGUE_VERSION"] = catalogue_version
    piezo["CATALOGUE_GRAMMAR"] = grammar
    piezo["PREDICTION_VALUES"] = values
    piezo["DRUG"] = drug
    piezo["SOURCE"] = json.dumps({})
    piezo["EVIDENCE"] = [
        json.dumps(
            {
                "solo_R": i[0][0],
                "solo_S": i[0][1],
                "background_R": i[1][0],
                "background_S": i[1][1],
                "p_value":i[2][0]
            }
        )
        if i
        else json.dumps({})
        for i in piezo["EVIDENCE"]
    ]
    piezo["OTHER"] = json.dumps({})

    piezo = piezo[
        [
            "GENBANK_REFERENCE",
            "CATALOGUE_NAME",
            "CATALOGUE_VERSION",
            "CATALOGUE_GRAMMAR",
            "PREDICTION_VALUES",
            "DRUG",
            "MUTATION",
            "PREDICTION",
            "SOURCE",
            "EVIDENCE",
            "OTHER",
        ]
    ]

    return piezo




def PiezoPredict(iso_df, catalogue_file, out_csv_matrix_file, out_text_error_file, out_text_summary, out_csv_predictions_file, drug, U_to_R=False, U_to_S=False, Print=True):
    catalogue = piezo.ResistanceCatalogue(catalogue_file)
    ids = iso_df.ena_run.unique().tolist()
    labels, predictions = [], []

    first_iteration = True  # Flag to track the first iteration

    for i in ids:
        df = iso_df[iso_df.ena_run == i]
        labels.append(df.phenotype.tolist()[0])
        mut_predictions = []
        for var in df.GENE_MUT:
            try:
                mut_predictions.append(catalogue.predict(var)[drug])
            except TypeError:
                mut_predictions.append("S")
            except ValueError:
                print(f"Something is up with {var}: No entry found in the catalogue for it")
                mut_predictions.append("S")
                mode = 'w' if first_iteration else 'a'  # Choose file mode based on first_iteration flag
                with open(out_text_error_file, mode) as error_file:
                    error_file.write(f"sample:{id}\tmut:{var}\tnote:Cannot find mut in catalogue\n")

        first_iteration = False  # Set the flag to False after the first iteration

        if "R" in mut_predictions:
            predictions.append("R")

        else:
            if "U" in mut_predictions:
                if U_to_R:
                    predictions.append("R")
                elif U_to_S:
                    predictions.append("S")
                else:
                    predictions.append("U")
            else:
                predictions.append("S")

    FN_id = []
    for i in range(len(labels)):
        if (predictions[i] == "S") & (labels[i] == "R"):
            FN_id.append(ids[i])

    # Create a dictionary from the data
    prediction_data = {'ids': ids, 'labels': labels, 'predictions': predictions}

    # Create a DataFrame
    prediction_df = pd.DataFrame(prediction_data)
    prediction_df.to_csv(out_csv_predictions_file)


    cm = confusion_matrix(labels, predictions)

    if "U" not in predictions:
        cm = cm[:2][:2]
    else:
        cm = cm[:2]
    # if print:
    #     print(cm)

    sensitivity = cm[0][0] / (cm[0][0] + cm[0][1])
    specificity = cm[1][1] / (cm[1][1] + cm[1][0])
    isolate_cov = (len(labels) - predictions.count("U")) / len(labels)

    if print:
        print("Catalogue coverage of isolates:", isolate_cov)
        print("Sensitivity:", sensitivity)
        print("Specificity:", specificity)

    # Open the output summary file in write mode
    with open(out_text_summary, 'w') as summary_file:
        # Write the summary content
        summary_file.write(f"This is a performance summary for {out_csv_predictions_file}.\n")
        summary_file.write(f"Catalogue coverage of isolates: {isolate_cov}\n")
        summary_file.write(f"Sensitivity: {sensitivity}\n")
        summary_file.write(f"Specificity: {specificity}\n")

    if predictions.count("U") == 0:
        df_cm = pd.DataFrame(cm, index=["R", "S"], columns=["R", "S"])
        df_cm.to_csv(out_csv_matrix_file)
        return [cm, isolate_cov, sensitivity, specificity, FN_id]

    if predictions.count("U") > 0:
        df_cm = pd.DataFrame(cm, index=["R", "S"], columns=["R", "S", "U"])
        df_cm.to_csv(out_csv_matrix_file)
        return [cm, isolate_cov, sensitivity, specificity, FN_id]




def DisplayPerformance(cm, columns):
    sns.set_context("talk")
    #cov = 100 * (len(val_samples) - cm[0][2] - cm[1][2]) / len(val_samples)
    df_cm = pd.DataFrame(cm, index=["R", "S"], columns=columns)
    plt.figure(figsize=(8, 5))
    sns.heatmap(
        df_cm, annot=True, cbar=False, fmt="g", cmap="Greens", annot_kws={"fontsize": 24}
    )
    plt.gca().invert_yaxis()
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.xlabel("Predicted", fontsize=14)
    plt.ylabel("True", fontsize=14)
    plt.show()

#%%
