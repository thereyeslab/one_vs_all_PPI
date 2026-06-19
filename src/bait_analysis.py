import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# read the data in tsv, create a df and filter only the interaction of the protein of interest, sort them by the requested column or columns and return the filtered df

def filter_proteins(aggregated_result_file, protein_of_interest, sort_columns):
    df_result = pd.read_csv(aggregated_result_file, sep="\t")
    df = df_result.copy()

    # filter the df to keep only the interactions of the protein of interest
    filtered_df = df[(df['protein_1'] == protein_of_interest) | (df['protein_2'] == protein_of_interest)]

    sorted_filtered_df = filtered_df.sort_values(by=sort_columns, ascending=False)

    return sorted_filtered_df


def plot_score_distribution(filtered_df, score_col):
    plt.figure(figsize=(8, 5))

    plt.hist(filtered_df[score_col].dropna(), bins=30)

    plt.xlabel(score_col)
    plt.ylabel("Number of interactions")
    plt.title(f"Distribution of {score_col}")

    plt.tight_layout()
    plt.savefig(f"score_distribution_{score_col}.png", dpi=300)
    plt.show()




def add_partner_column(df, protein_of_interest):
    """
    Add a column containing the interaction partner of the protein of interest.
    """
    df = df.copy()
    df["partner"] = df.apply(
        lambda row: row["protein_2"] if row["protein_1"] == protein_of_interest else row["protein_1"],
        axis=1
    )
    return df

def plot_ranked_score_dots(filtered_df, protein_of_interest, score_col, top_n=50):
    df = add_partner_column(filtered_df, protein_of_interest)
    df = df.sort_values(score_col, ascending=False).head(top_n)

    plt.figure(figsize=(10, 10))
    plt.scatter(df[score_col], df["partner"])

    plt.gca().invert_yaxis()
    plt.xlabel(score_col)
    plt.ylabel("Interaction partner")
    plt.title(f"Ranked interaction scores for {protein_of_interest}")

    plt.tight_layout()
    plt.savefig(f"ranked_interaction_scores_{protein_of_interest}.png", dpi=300)
    plt.show()
    # save the plot as a png file
    

if __name__ == "__main__":
    aggregated_result_file = "/Volumes/Masoumeh/af3/af3_output_zip/pair_scores_size_corrected_aggregated.tsv"
    protein_of_interest = "CCNA_01764"  # replace with your protein of interest
    sort_columns = ["raw_score_resolved"]  # replace with the columns you want to sort by

    filtered_df = filter_proteins(aggregated_result_file, protein_of_interest, sort_columns)
    # save the filtered df to a new tsv file
    output_file = f"filtered_{protein_of_interest}_{'_'.join(sort_columns)}.tsv"
    filtered_df.to_csv(output_file, sep="\t", index=False)
    print(f"Filtered results saved to {output_file}")


    plot_ranked_score_dots(
    filtered_df,
    protein_of_interest="CCNA_01764",
    score_col="raw_score_resolved",
    top_n=50
    )



    plot_score_distribution(filtered_df, "raw_score_resolved")


