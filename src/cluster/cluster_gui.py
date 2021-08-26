from bokeh.layouts import grid, layout, column
from bokeh.models import CustomJS, ColumnDataSource, Button, HoverTool
from bokeh.models.widgets import DataTable, TableColumn, Slider
from bokeh.client import push_session
from bokeh.tile_providers import CARTODBPOSITRON, get_provider
from bokeh.models import LinearColorMapper, TextInput
import numpy as num
from pyrocko import model
from bokeh.plotting import figure, output_file, curdoc, show
import math
import copy
import sys
from bokeh.server.server import Server
from matplotlib.pyplot import cm
from bokeh.palettes import Dark2_5 as palette
# itertools handles the cycling
import itertools

def list_duplicates_of(seq, item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item, start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs


def button_callback():
    server.io_loop.stop()
    server.stop()  # Stop the server


def clearbutton_callback():
    s2.data = {k: [] for k in s2.data}


def clearbuttonlast_callback():
    last_selection = len(s1.selected.indices)
    last_table_len = len(s2.data["x"])
    s2.data = {"x": s2.data["x"][0:-last_selection], "y": s2.data["y"][0:-last_selection], "ind": s2.data["ind"][0:-last_selection], "name": s2.data["name"][0:-last_selection]}

def clearbuttonone_callback():
    last_selection = 1
    s2.data = {"x": s2.data["x"][0:-last_selection], "y": s2.data["y"][0:-last_selection], "ind": s2.data["ind"][0:-last_selection], "name": s2.data["name"][0:-last_selection]}


def update():
    global ds, old_indices
    old_indices_current = copy.deepcopy(old_indices)
    print(slider.value)
    print(text_input.value)
    slider.value = int(text_input.value)
    #slider.value = int(text_input.value)
#    for name in s2.data["name"]:
#        idx_duplicates = list_duplicates_of(s2.data["name"], name)
    if len(old_indices_current) != len(s1.selected.indices):
        old_indices = s1.selected.indices
        ds = ColumnDataSource(data=dict(x=[], y=[], ind=[], name=[]))
        x_stream = []
        y_stream = []
        name_stream = []
        ind_stream = []
        for ind in s1.selected.indices:
            x_stream.append(s1.data["x"][ind])
            y_stream.append(s1.data["y"][ind])
            name_stream.append(s1.data["name"][ind])
            ind_stream.append(s1.data["ind"][ind])
        ds.stream({"x": x_stream, "y": y_stream, "ind": ind_stream, "name": name_stream})
        tr_data = []
        tr_xdata = []
        colors = itertools.cycle(palette)
        s_tr_ds = []
        compare_ds = copy.deepcopy(ds)
        plotted = []
        fig03.renderers = []
        for i, tr in enumerate(traces):
            for name in compare_ds.data["name"]:
                if tr.network+"."+tr.station+"." == name:
                    if name not in plotted:
                        plotted.append(name)
                        name_vector = num.repeat(name, len(tr.ydata))
                        ind_vector = num.repeat(i, len(tr.ydata))

                        s_tr = ColumnDataSource(data=dict(x=tr.ydata, y=event.time-tr.get_xdata(), name=name_vector, ind=ind_vector))
                        fig03.line("y", "x", source=s_tr, alpha=1, color=next(colors))
    cluster_result = copy.deepcopy(s2.data)


def merc(lat, lon):
    r_major = 6378137.000
    x_m = r_major * math.radians(lon)
    scale = x_m/lon
    y_m = 180.0/math.pi * math.log(math.tan(math.pi/4.0 +
                                   lat * (math.pi/180.0)/2.0)) * scale
    return (x_m, y_m)


def cluster_gui(doc):
    global s2, s1, old_indices
    old_indices = []
    output_file("tile.html")
    tile_provider = get_provider(CARTODBPOSITRON)
    x = []
    y = []
    name = []
    global fig03
    fig03 = figure(
        plot_width=400,
        plot_height=400,
        tools=["box_zoom", "wheel_zoom", "reset", "save"],
        title="Waveforms from current selection",
    )
    for i, st in enumerate(stations):
        xm, ym = merc(st.lat, st.lon)
        x.append(xm)
        y.append(ym)
        name.append(st.nsl_string())


    # create first subplot
    plot_width = 400
    plot_height = 400
    d = num.ones_like(x)

    s1 = ColumnDataSource(data=dict(x=x, y=y, ind=d, name=name))
    # range bounds supplied in web mercator coordinates
    fig01 = figure(x_axis_type="mercator", y_axis_type="mercator",
                   plot_width=plot_width,
                   plot_height=plot_height,
                   tools=["lasso_select", "box_select", "reset",
                          "save", "box_zoom", "wheel_zoom"],
                   title="Select",)
    fig01.add_tile(tile_provider)


    fig01.scatter("x", "y", source=s1, alpha=0.6, size=8)

    # create second subplot
    s2 = ColumnDataSource(data=dict(x=[], y=[], ind=[], name=[]))

    color_mapper = LinearColorMapper(palette='Magma256', low=1, high=100)

    fig02 = figure(
        x_axis_type="mercator", y_axis_type="mercator",
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=(num.min(x), num.max(x)),
        y_range=(num.min(y), num.max(y)),
        tools=["box_zoom", "wheel_zoom", "reset", "save"],
        title="Stations selected for Array",
    )
    fig02.add_tile(tile_provider)


    fig02.scatter("x", "y", source=s2, alpha=1,
                  color={'field': 'ind', 'transform': color_mapper}, size=8)

    x_event, y_event = merc(event.lat, event.lon)
    fig01.scatter(x_event, y_event, size=8, color="red")
    fig02.scatter(x_event, y_event, size=8, color="red")

    columns = [
        TableColumn(field="x", title="X axis"),
        TableColumn(field="y", title="Y axis"),
        TableColumn(field="ind", title="indices"),
        TableColumn(field="name", title="name"),
    ]

    table = DataTable(
        source=s2,
        columns=columns,
        width=400,
        height=600,
        sortable=True,
        selectable=True,
        editable=True,

    )

    source_count = 0
    callback_slider = CustomJS(code="""
        source_count = slider.value;
        """)
    global slider

    slider = Slider(start=1, end=100, value=1, step=1, title="Array number")
    slider.js_on_change('value', callback_slider)

    s1.selected.js_on_change(
        "indices",
        CustomJS(
            args=dict(s1=s1, s2=s2, s3=slider, table=table),
            code="""
            var inds = cb_obj.indices;
            var d1 = s1.data;
            var d2 = s2.data;
            const A = s3.value;

            for (var i = 0; i < inds.length; i++) {
                d2['x'].push(d1['x'][inds[i]])
                d2['y'].push(d1['y'][inds[i]])
                d2['name'].push(d1['name'][inds[i]])
                d2['ind'].push(A)
            }
            s2.change.emit();
            table.change.emit();
    	    s2.data = s2.data;

            var inds = source_data.selected.indices;
            var data = source_data.data;
            var out = "name, x, y, ind\\n";
            for (i = 0; i < inds.length; i++) {
                out += data['name'][inds[i]] + "," + data['x'][inds[i]] + "," + data['y'][inds[i]] + "," + data['ind'][inds[i]]   + "\\n";
            }
            var file = new Blob([out], {type: 'text/plain'});

        """),
        ),


    savebutton = Button(label="Save", button_type="success")
    savebutton.callback = CustomJS(
        args=dict(source_data=s1),
        code="""
            var inds = source_data.selected.indices;
            var data = source_data.data;
            var out = "name, x, y, ind\\n";
            for (i = 0; i < inds.length; i++) {
                out += data['name'][inds[i]] + "," + data['x'][inds[i]] + "," + data['y'][inds[i]] + "," + data['ind'][inds[i]]   + "\\n";
            }
            var file = new Blob([out], {type: 'text/plain'});
            var elem = window.document.createElement('a');
            elem.href = window.URL.createObjectURL(file);
            elem.download = 'arrays.txt';
            document.body.appendChild(elem);
            elem.click();
            document.body.removeChild(elem);
            """,
    )


    tooltips = [
        ("X:", "@x"),
        ("Y:", "@y"),
        ("Array:", "@ind"),
        ("Station:", "@name"),

    ]

    fig01.add_tools(HoverTool(tooltips=tooltips))
    fig02.add_tools(HoverTool(tooltips=tooltips))
    fig03.add_tools(HoverTool(tooltips=tooltips))

    endbutton = Button(label="End and proceed", button_type="success")
    endbutton.on_click(button_callback)

    clearbutton = Button(label="Clear all", button_type="success")
    clearbutton.on_click(clearbuttonlast_callback)
    clearbuttonlast = Button(label="Clear last selection", button_type="success")
    clearbuttonlast.on_click(clearbuttonlast_callback)
    clearbuttonone = Button(label="Remove one from list", button_type="success")
    clearbuttonone.on_click(clearbuttonone_callback)

    b = Button(label="Reset all plots")
    b.js_on_click(CustomJS(code="""\
document.querySelectorAll('.bk-tool-icon-reset[title="Reset"]').forEach(d => d.click())
"""))
    #layout = grid([fig01, fig02, table, fig03, slider, clearbuttonlast, clearbutton, savebutton, endbutton], ncols=3, nrows=4)
    global text_input
    text_input = TextInput(value="1", title="Array number:")

    buttons = column(clearbuttonlast, clearbuttonone, clearbutton, savebutton, endbutton)
    inputs = column(text_input, slider)

    layout_grid = layout([fig01, fig02, buttons],[fig03, inputs, table])

    #curdoc().add_root(layout)
    cluster_result = []
    doc.add_root(layout_grid)
    #global session
    #session = push_session(curdoc())

    #curdoc().add_periodic_callback(update, 100)
    #session.show(layout)
    doc.add_periodic_callback(update, 900)
    curdoc().title = "Array selection"
    #session.loop_until_closed()
    #show(layout)

    #show(layout)


def start_server(stations_in, traces_in, event_in):
    global stations
    stations = stations_in
    global traces
    traces = traces_in
    global event
    event = event_in
    global server
    server = Server({'/': cluster_gui}, num_procs=1)
    server.start()
    server.io_loop.add_callback(server.show, "/", )
    server.io_loop.start()

    return ds

