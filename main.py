from dash import Dash, dcc, html, Input, Output
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go

import helper

app = Dash(__name__)

#################
# Default layout settings:
global_layout = go.Layout(
    template="plotly_white",
    font={"family": "sans-serif",
          "size": 15,
          "color": "black"},
    xaxis={"showgrid": False,
           "showline": True,
           "linecolor": "black",
           "ticks": "inside",
           "mirror": True},    
    yaxis={"showgrid": False,
           "showline": True,
           "linecolor": "black",
           "ticks": "inside",
           "mirror": True},  
    showlegend=True,
    legend={"x": 0,
            "y": 1,
            "traceorder": "normal",
            "font": {"family": "sans-serif",
                     "size": 12,
                     "color": "black"}},
)

common_style_params = {"display": "inline-block",
                       "vertical-align": "top",
                       "text-align": "center",
                       "margin": "0",
                       "padding": "0"}

text_style_params = {"text-align": "center",
                     "color": "#85c1e9"}
#################
# Load data
GENOME =  pd.read_csv("data/genes.csv", index_col=0).T.to_dict()
GENOME = {gene: (val["0"], val["1"], val["2"]) for gene, val in GENOME.items()}
GENOME.update({"ori": (0, 10, "none")})

DF_PROTEOMICS = pd.read_csv("data/schmidt2016_s6.csv")

MAP_GENE_PROTEIN = DF_PROTEOMICS[["Gene", "Uniprot Accession"]].set_index("Gene").to_dict()["Uniprot Accession"]


##################
# App layout 
app.layout = html.Div(
    style={"text-align": "center"},
    children=[ 
    html.H1([html.I("E. coli"), " BW25113 genome and proteome"], style=text_style_params),

    html.Br(),

    html.Div([
        html.Div([html.H2("Gene location", style=text_style_params),
                  dcc.Dropdown(list(MAP_GENE_PROTEIN.keys()), value="aas", id="gene_selector", multi=True, style={"width": "95%"}),
                  dcc.Graph(id="chromosome", figure={}, style={"width": "95%"}),
                  ], 
                style = {"width": "32%", **common_style_params}
                ),
        
        html.Div([html.H2("Protein abundance across conditions", style=text_style_params),
                  dcc.Dropdown(["Absolute", "Normalized"], value="Absolute", id="protein_values", style={"width": "95%"}),
                  dcc.Graph(id="proteomics", figure={}, style={"width": "95%"}),
                  ], 
                 style = {"width": "58%", **common_style_params},
                 ),                 
        ],
        style={"margin": "0", "padding": "0"}
        ),

    html.Div([
        html.Div([html.H2("Heatmaps", style=text_style_params),
                  dcc.Dropdown(["gene-gene distance", "proteome correlation"], value="proteome correlation", id="choice_graph", style={"width": "95%"}),
                  dcc.Graph(id="choice", figure={}),
                  ], 
                style = {"width": "29%", **common_style_params}
                ),
        
        html.Div([html.H2("Information", style=text_style_params),
                  html.Div("(click on a cell from the heatmap)", style={"text-align": "center", "vertical-align": "top"}),
                  html.Iframe(id="website1", src="", style={"width": "49%", "height": "600px", "border": "2px black solid"}),
                  html.Iframe(id="website2", src="", style={"width": "49%", "height": "600px", "border": "2px black solid"}),
                  ], 
                 style = {"width": "65%", **common_style_params},
                 ),
        ],
        style={"margin": "0", "padding": "0"}
        ),
],  
)



#################
# Callbacks
@app.callback(
    Output("chromosome", "figure"),
    Output("proteomics", "figure"),
    Output("choice", "figure"),
    [Input("gene_selector", "value"),
     Input("protein_values", "value"),
     Input("choice_graph", "value")]
)
def update_chromosome(value, prot_values, choice):

    palette = px.colors.qualitative.Plotly

    # Display chromosome:
    fig = helper.plot_circle(value)

    to_display = ["ori"] 
    if value is not None:
        to_display += value

    # Find coordinates of genes:
    colors_gene = {}
    i = 0
    for gene in to_display:

        pos_start, pos_end, orientation = GENOME[gene]
        theta1, theta2 = helper.find_coordinates_gene(pos_start, pos_end, genome_length=4.6e6)

        if gene in ["ori"]:
            c = "yellow"
        else:
            c = palette[i]
            i += 1

        fig = helper.add_arc(fig, theta1, theta2, gene, orientation, color=c) 
        colors_gene[gene] = c

    # Select proteomics data:
    proteins_to_display = [MAP_GENE_PROTEIN[gene] for gene in value]
    colors_prot = {MAP_GENE_PROTEIN[gene]: colors_gene[gene] for gene in value}

    df_plot = helper.filter_proteomics(DF_PROTEOMICS, proteins_to_display, prot_values)
    fig2 = helper.display_proteomics(df_plot, colors=colors_prot)
    if prot_values == "Absolute":
        fig2.update_layout(yaxis_title="Protein abundance [copies/cell]",
                           xaxis_title="Condition")
    elif prot_values == "Normalized":
        fig2.update_layout(yaxis_title="Normalized protein abundance [-]",
                           xaxis_title="Condition")

    # Normalize proteomics data across conditions:
    df_pivot = df_plot.pivot(columns="Uniprot Accession", index="condition", values="ProteinAbundance")
    if prot_values == "Absolute":
        df_pivot = df_pivot.transform(lambda x: (x - x.mean()) / x.std(), axis=0)

    # Correlations between proteins:
    corr = df_pivot.corr()

    # Gene-to-gene distance:
    distance_matrix = helper.calc_gene_distance(value, GENOME, genome_length=4.6e6)
    distance_matrix = helper.make_symmetric(distance_matrix.to_numpy())

    if choice == "proteome correlation":  # correlation between proteome trends
        fig3 = px.imshow(corr.to_numpy(),
                        x=list(corr.index), y=list(corr.index),
                        title="Correlation between protein abundances",
                        color_continuous_scale="RdBu", zmin=-1, zmax=1,
                        text_auto=".2f")
        fig3 = helper.update_label_colors(fig3, list(corr.index), colors_prot)

    elif choice == "gene-gene distance":  # gene-to-gene distance
        fig3 = px.imshow(distance_matrix, x=value, y=value,
                        title="Gene-to-gene distance",
                        color_continuous_scale="Blues", zmin=0,
                        text_auto=True)
        fig3 = helper.update_label_colors(fig3, value, colors_gene)


    fig.update_layout(global_layout)
    fig2.update_layout(global_layout)
    fig3.update_layout(global_layout)

    return fig, fig2, fig3

@app.callback(
    Output("website1", "src"),
    Output("website2", "src"),
    [Input("choice", "clickData"),
     Input("choice_graph", "value")])
def retrieve_websites(clickData, choice):

    base_url = "https://www.uniprot.org/uniprotkb/"
    blank = "about:blank"
    
    if clickData is None:
        return blank, blank

    clicked_point = clickData["points"][0]

    if choice == "proteome correlation":
        prot1 = clicked_point["x"]
        prot2 = clicked_point["y"]
    elif choice == "gene-gene distance":
        gene1 = clicked_point["x"]
        gene2 = clicked_point["y"]
        prot1 = MAP_GENE_PROTEIN[gene1]
        prot2 = MAP_GENE_PROTEIN[gene2]
    else:
        prot1 = ""
        prot2 = ""


    if prot1 != prot2:
        return base_url+prot1, base_url+prot2
    else:
        return base_url+prot1, blank


if __name__ == "__main__":
    app.run(debug=True)