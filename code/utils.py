import csv
import pandas as pd
from sklearn.metrics import confusion_matrix
import piezo
import matplotlib.pyplot as plt
import seaborn as sns
import sys

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
    if not set(df.columns.str.upper()) == required_columns:
        raise ValueError("Genes file must have the following columns: drug, gene")


def check_mutations_file_columns(df):
    """
   Check that the columns in the dataframe from mutations file are correct
   """
    required_columns = {"UNIQUEID", "GENE", "MUTATION"}
    if not set(df.columns.str.upper()) == required_columns:
        raise ValueError("Mutations file must have the following columns: UNIQUEID, GENE, MUTATION")

def check_phenotypes_file_columns(df):
    """
   Check that the columns in the dataframe from phenotypes file are correct
   """
    required_columns = {"UNIQUEID", "DRUG", "PHENOTYPE", "PHENOTYPE_QUALITY"}
    if not set(df.columns.str.upper()) == required_columns:
        raise ValueError("Phenotypes file must have the following columns: UNIQUEID, DRUG, PHENOTYPE, PHENOTYPE_QUALITY")


def get_genes_of_interest_from_genes_file(filename, drug):
    """
    Get list of genes by filtering genes list df by drug of interest
    """
    genes_separator = get_separator(filename)
    genes_df = pd.read_csv(filename, sep=genes_separator)
    check_gene_file_columns(genes_df)

    genes_gene_column = next((col for col in genes_df.columns if 'gene' in col.lower()), None)
    genes_drug_column = next((col for col in genes_df.columns if 'drug' in col.lower()), None)

    if genes_gene_column and genes_drug_column:
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
        if "R" in group["PHENOTYPE"].values:
            return group[group["PHENOTYPE"] == "R"].iloc[0:1]
        else:
            return group.iloc[0:1]
    else:
        return group

def RSIsolateTable(df, genes):
    """returns df of number of R and S isolates"""
    table = {}
    table["Total"] = {
        "R": df[df.PHENOTYPE == "R"].UNIQUEID.nunique(),
        "S": df[df.PHENOTYPE == "S"].UNIQUEID.nunique(),
        "Total": df.UNIQUEID.nunique(),
    }
    for i in genes:
        d = df[df.GENE == i]
        table[i] = {
            "R": d[d.PHENOTYPE == "R"].UNIQUEID.nunique(),
            "S": d[d.PHENOTYPE == "S"].UNIQUEID.nunique(),
            "Total": d[d.PHENOTYPE == "R"].UNIQUEID.nunique()
            + d[d.PHENOTYPE == "S"].UNIQUEID.nunique(),
        }

    return pd.DataFrame.from_dict(table).T


def RSVariantTable(df, genes):
    """returns df of number of R and S variants"""
    table = {}
    table["Total"] = {
        "R": df[df.PHENOTYPE == "R"].UNIQUEID.count(),
        "S": df[df.PHENOTYPE == "S"].UNIQUEID.count(),
        "Total": df.UNIQUEID.count(),
    }
    for i in genes:
        d = df[df.GENE == i]
        table[i] = {
            "R": d[d.PHENOTYPE == "R"].UNIQUEID.count(),
            "S": d[d.PHENOTYPE == "S"].UNIQUEID.count(),
            "Total": d[d.PHENOTYPE == "R"].UNIQUEID.count()
            + d[d.PHENOTYPE == "S"].UNIQUEID.count(),
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



def PrepTables(drug, genes, train_ids, valid_ids):
    mutations = pd.read_csv(f"data/demo/{drug.lower()}_mutations.csv")
    genomes = pd.read_pickle("data/cryptic-tables-v2.0.1/GENOMES.pkl.gz").reset_index()

    phenotypes = pd.read_pickle(
        "data/cryptic-tables-v2.0.1/DST_MEASUREMENTS.pkl.gz"
    ).reset_index()
    phenotypes = phenotypes[
        (phenotypes.DRUG == drug) & (phenotypes.PHENOTYPE.isin(["R", "S"]))
    ]
    # if multiple phenos, keep R else first
    phenotypes = (
        phenotypes.groupby("UNIQUEID")
        .apply(filter_multiple_phenos)
        .reset_index(drop=True)
    )

    samples = pd.merge(genomes, phenotypes, on=["UNIQUEID"], how="inner")
    samples = samples[["UNIQUEID", "DRUG", "PHENOTYPE", "METHOD_MIC", "QUALITY"]]

    train_samples = samples[
        samples.UNIQUEID.isin(train_ids) & (samples.QUALITY.isin(["HIGH", "MEDIUM"]))
    ]
    val_samples = samples[samples.UNIQUEID.isin(valid_ids)]

    mutations = mutations[mutations.GENE.isin(genes)]

    mutations["GENE_MUT"] = [
        f"{row['GENE']}@{row['MINOR_MUTATION'] if row['IS_MINOR_ALLELE'] else row['MUTATION']}"
        for _, row in mutations.iterrows()
    ]

    mutations["IS_SYNONYMOUS"] = [
        row["MUTATION"][0] == row["MUTATION"][-1] for _, row in mutations.iterrows()
    ]

    mutations["FRS"] = [
        1 if ~mutations["IS_MINOR_ALLELE"][i] else mutations["FRS"][i]
        for i in mutations.index
    ]

    train_mutations = mutations[
        (mutations.UNIQUEID.isin(train_ids)) & (~mutations.IS_SYNONYMOUS)
    ]
    val_mutations = mutations[(mutations.UNIQUEID.isin(valid_ids))]

    return train_samples, val_samples, train_mutations, val_mutations


def PiezoPredict(iso_df, catalogue_file, drug, U_to_R=False, U_to_S=False, Print=True):
    catalogue = piezo.ResistanceCatalogue(catalogue_file)

    ids = iso_df.UNIQUEID.unique().tolist()
    labels, predictions = [], []

    for i in ids:
        df = iso_df[iso_df.UNIQUEID == i]
        labels.append(df.PHENOTYPE.tolist()[0])
        mut_predictions = []
        for var in df.GENE_MUT:
            try:
                mut_predictions.append(catalogue.predict(var)[drug])
            except TypeError:
                mut_predictions.append("S")

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

    cm = confusion_matrix(labels, predictions)

    if "U" not in predictions:
        cm = cm[:2][:2]
    else:
        cm = cm[:2]
    if print:
        print(cm)

    sensitivity = cm[0][0] / (cm[0][0] + cm[0][1])
    specificity = cm[1][1] / (cm[1][1] + cm[1][0])
    isolate_cov = (len(labels) - predictions.count("U")) / len(labels)

    if print:
        print("Catalogue coverage of isolates:", isolate_cov)
        print("Sensitivity:", sensitivity)
        print("Specificity:", specificity)

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
