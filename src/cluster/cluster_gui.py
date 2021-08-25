from bokeh.layouts import grid
from bokeh.models import CustomJS, ColumnDataSource, Button, HoverTool
from bokeh.models.widgets import DataTable, TableColumn, Slider
from bokeh.client import push_session
from bokeh.tile_providers import CARTODBPOSITRON, get_provider
from bokeh.models import LinearColorMapper

import numpy as num
from pyrocko import model
from bokeh.plotting import figure, output_file, curdoc
import math

output_file("tile.html")
tile_provider = get_provider(CARTODBPOSITRON)


def merc(lat, lon):
    r_major = 6378137.000
    x_m = r_major * math.radians(lon)
    scale = x/lon
    y_m = 180.0/math.pi * math.log(math.tan(math.pi/4.0 +
                                 lat * (math.pi/180.0)/2.0)) * scale
    return (x_m, y_m)


stations = model.load_stations("stations.txt")
x = []
y = []
name = []
for i, st in enumerate(stations):
    xm, ym = merc(st.lat, st.lon)
    x.append(xm)
    y.append(ym)
    name.append(st.nsl_string())

# create first subplot
plot_width = 600
plot_height = 600
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
count = 0

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

slider = Slider(start=1, end=100, value=1, step=1, title="Array")
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
]

fig02.add_tools(HoverTool(tooltips=tooltips))

layout = grid([fig01, fig02, table, savebutton, slider], ncols=3, nrows=2)
cluster_result = []

def update():
    cluster_result = s2.data
    print(s2.data)

session = push_session(curdoc())
curdoc().add_periodic_callback(update, 100)
session.show(layout)
session.loop_until_closed()

#show(layout)
