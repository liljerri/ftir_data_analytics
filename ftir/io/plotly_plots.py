from plotly.graph_objs import Layout
import plotly.io as pio
import plotly.graph_objects as go

# Method to define plotly.graph_objects plotting format


class plotly_go:
    def __init__(
        self, left=False, right=False, min_y=0, max_y=100
    ):  # Runs every time we run plotly_go

        self.left = left  # Self refers to the instance, i.e. emp_1.first
        self.right = right
        self.min_y = min_y
        self.max_y = max_y
        self.layout = Layout(
            paper_bgcolor="white",
            plot_bgcolor="white",
            xaxis={
                "showgrid": False,
                "showline": True,
                "linewidth": 2,
                "linecolor": "black",
                "mirror": True,
                "title_text": "Wavenumber",
                "range": [self.left, self.right],
                "ticks": "outside",
                "tickcolor": "black",
                "tickwidth": 2,
                "ticklen": 10,
            },
            yaxis={
                "showgrid": False,
                "showline": True,
                "linewidth": 2,
                "linecolor": "black",
                "mirror": True,
                "title_text": "Abs",
                # "range": [self.min_y, self.max_y],
                "ticks": "outside",
                "tickcolor": "black",
                "tickwidth": 2,
                "ticklen": 10,
            },
        )


class raw_data_graph(plotly_go):
    def __init__(
        self, df
    ):  # Never pass mutable args!
        super().__init__()  # Allows Employee class to pass info
        
        

        # Create interactive graph of cgms

        # Call graph object figure initialization
        fig = go.Figure(layout=self.layout)
        fig.update_layout(xaxis=dict(tickmode="linear", tick0=0, dtick=2))
        # Add traces
        for i in range(len(df.columns) - 1):

            fig.add_trace(
                go.Scatter(
                    x=df.iloc[:, 0],
                    y=df.iloc[:, (i + 1)],
                    mode="lines",
                    name=df.columns[(i + 1)],
                )
            )

        pio.show(fig)

class baseline_check_graph(plotly_go):
    def __init__(
        self, x_val, y_val, base,left=False, right=False, min_y=0, max_y=100
    ):
        super().__init__(
            left=False, right=False, min_y=0, max_y=100
        )  # Allows Employee class to pass info
        self.base = base
        # Call graph object figure initialization
        fig = go.Figure(layout=self.layout)

        # Add traces
        fig.add_trace(go.Scatter(x = x_val, y=(y_val - self.base),
                            mode='lines',
                            name='baseline corrected ref1'))
        fig.add_trace(go.Scatter(x = x_val, y=y_val-8.7,
                            mode='lines',
                            name='ref1'))
        fig.add_trace(go.Scatter(x = x_val, y=self.base-8.7,
                            mode='lines',
                            name='baseline'))

        pio.show(fig)    

class peak_ctr_check_graph(plotly_go):
    def __init__(
        self, x_val, y_val, indexes, left=False, right=False, min_y=0, max_y=100
    ):
        super().__init__(
            left=False, right=False, min_y=0, max_y=100
        )  # Allows Employee class to pass info
        self.indexes = indexes
        # Call graph object figure initialization
        fig = go.Figure(layout=self.layout)

        # Add traces
        fig.add_trace(go.Scatter(x = x_val, y=(y_val),
                            mode='lines',
                            name='baseline corrected ref1'))
        fig.add_trace(go.Scatter(x = x_val[self.indexes], y= y_val[self.indexes],
                            mode='markers',
                            marker={'size':8, 'color':'#EB55BF'}))
                                
        pio.show(fig) 
