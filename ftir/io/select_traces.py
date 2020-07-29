'''from: https://gist.github.com/pbugnion/5bb7878ff212a0116f0f1fbc9f431a5c and https://stackoverflow.com/questions/57219796/ipywidgets-dynamic-creation-of-checkboxes-and-selection-of-data'''
import difflib
import random
import requests
import ipywidgets as widgets

from ftir.io.plotly_plots import raw_data_graph

import ipywidgets as widgets
import itertools

def multi_checkbox_widget(sample_names):
    flatten = itertools.chain.from_iterable #tool to flatten a nested list
    names = []
    checkbox_objects = []
    label_objects = []
    for name in sample_names:
        checkbox_objects.append(widgets.Checkbox(layout = widgets.Layout(width='50%')))
        label_objects.append(widgets.Label(name,layout = widgets.Layout(width='25%')))
        names.append(name)
    arg_dict = {names[i]: checkbox for i, checkbox in enumerate(checkbox_objects)}
    labels_checkboxes =  list(flatten(zip(label_objects, checkbox_objects)))
    ui = widgets.HBox(labels_checkboxes, layout=widgets.Layout(flex_flow='row wrap'))
    selected_data = []
    def select_data(**kwargs):
        selected_data.clear()
        for key in kwargs:
            if kwargs[key] is True:
                selected_data.append(key)
    display(ui)
    widgets.interactive_output(select_data, arg_dict)
    return(selected_data)

def selected_files_df(selection, df):
    selection.insert(0, df.columns[0])
    raw_data_graph(df[selection])
    return df[selection].copy()



def multi_checkbox_grapher(sample_names, df):
    import ipywidgets as widgets
    import itertools
    from ftir.io.plotly_plots import raw_data_graph
    output = widgets.Output()
    plot_output = widgets.Output()
    
    flatten = itertools.chain.from_iterable #tool to flatten a nested list
    names = []
    checkbox_objects = []
    label_objects = []
    for name in sample_names:
        checkbox_objects.append(widgets.Checkbox(layout = widgets.Layout(width='50%')))
        label_objects.append(widgets.Label(name,layout = widgets.Layout(width='25%')))
        names.append(name)
    arg_dict = {names[i]: checkbox for i, checkbox in enumerate(checkbox_objects)}
    labels_checkboxes =  list(flatten(zip(label_objects, checkbox_objects)))
    ui = widgets.HBox(labels_checkboxes, layout=widgets.Layout(flex_flow='row wrap'))
    selected_data = []
    def select_data(**kwargs):
        plot_output.clear_output()
        with output:
            selected_data.clear()
            for key in kwargs:
                if kwargs[key] is True:
                    selected_data.append(key)
        with plot_output:
            selected_data.insert(0, df.columns[0])
            raw_data_graph(df[selected_data])

    display(ui)
    widgets.interactive_output(select_data, arg_dict)
    display(plot_output)
    return(selected_data)
