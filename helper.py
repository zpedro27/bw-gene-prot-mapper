import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
import itertools

def plot_circle(title=""):
    
    x, y = _get_arc_coordinates(0, 2*np.pi, r=1, npoints=100)

    # Create the figure
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=x, y=y,
        name="chromosome",
        showlegend=False,
        line=dict(color="lightgray", width=10)
    ))
    fig.update_layout(
        xaxis=dict(scaleanchor="y", visible=False),
        yaxis=dict(visible=False),
        width=450,
        height=450
    )
    return fig

def add_arc(fig, theta1, theta2, gene, orientation, color):
    x, y = _get_arc_coordinates(theta1, theta2, r=1, npoints=10)

    fig.add_trace(go.Scatter(
        x=x, y=y,
        name=gene,
        line=dict(width=25, color=color),
        opacity=0.4,
    ))
    return fig

def _get_arc_coordinates(theta1, theta2, r=1, npoints=100):
    theta = np.linspace(theta1, theta2, npoints)
    x = r * np.cos(theta)
    y = r* np.sin(theta)
    return x, y

def find_coordinates_gene(pos_start, pos_end, genome_length):
    theta1 = 2 * np.pi * pos_start / genome_length
    theta2 = 2 * np.pi * pos_end / genome_length
    return theta1, theta2

def display_proteomics(df, colors: dict):
    fig = px.scatter(df, x="condition", y="ProteinAbundance",
             color="Uniprot Accession",
            color_discrete_map=colors,
             height=400)
    for x in fig.data:
        x.update(mode="markers+lines")
    return fig

def filter_proteomics(df, proteins, value_type):

    df_out = df.loc[df["Uniprot Accession"].isin(proteins)].sort_values(by="condition")

    if value_type == "Normalized":
        df_out["ProteinAbundance"] = df_out.groupby(["Uniprot Accession"])["ProteinAbundance"].transform(lambda x: (x - x.mean()) / x.std())
        return df_out
    elif value_type == "Absolute":
        return df_out
    else:
        raise ValueError(f"Unknown value_type={value_type}")
    return None

def calc_gene_distance(gene_list, genome, genome_length):

    all_pairs = itertools.combinations_with_replacement(gene_list, 2)

    distances = {}
    for pair in all_pairs:
        gene1, gene2 = pair
        coords1 = genome[gene1]
        coords2 = genome[gene2]
        distances[pair] = _calc_gene_distance(coords1[0:2], coords2[0:2])
    return _convert_paired_distances_to_df(distances)

def _convert_paired_distances_to_df(distances: dict):
    df = pd.DataFrame(distances, index=[0]).T
    df = df.unstack()
    df.columns = df.columns.droplevel(0)
    return df

def _calc_gene_distance(coords1: tuple, coords2: tuple):
    all_combs = itertools.product(coords1, coords2)
    all_angular_distances = [abs(x[0] - x[1]) for x in all_combs]
    return min(all_angular_distances)


def make_symmetric(array: np.array):
    return np.triu(array) + np.triu(array, 1).T


def update_label_colors(fig, labels, colors):

    fig.update_xaxes(
        tickvals=[i for i in range(len(labels))],
        ticktext=[
            f'<span style="color:{colors[label]};">{label}</span>'
            for label in labels
        ]
    )

    fig.update_yaxes(
        tickvals=[i for i in range(len(labels))],
        ticktext=[
            f'<span style="color:{colors[label]};">{label}</span>'
            for label in labels
        ]
    )

    return fig