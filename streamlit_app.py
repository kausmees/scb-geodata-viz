import streamlit as st
import pandas as pd
import geopandas as gpd
import pydeck as pdk
import numpy as np
import math
from shapely import geometry

st.set_page_config(page_title="Statistics Sweden open geodata", page_icon=":lemon:", layout="wide")

def transform_df_coords(df, transf):
    """
    Transform the coordinates of a geopandas dataframe containing polygons.

    @param df: the dataframe
    @param transf: function to transform coordinates with
    @return: new dataframe with transformed coordinates
    """
    for index, row in df.iterrows():
        xx, yy = row.geometry.exterior.coords.xy
        new_coords = list(map(transf, xx ,yy))
        new_polygon = geometry.Polygon(new_coords)
        row.geometry = new_polygon
        df.iloc[index] = row
    return df

def grid_to_geodetic(y, x) :
    """
    Conversion from SWEREF99TM to lat, lon.
    Using code from here https://gis.stackexchange.com/questions/292637/from-sweref99-to-epsg3857

    @param y: North coordinate
    @param x: East coordinate
    @return: [lat, lon]
    """
    def math_sinh(value):
        return 0.5 * (math.exp(value) - math.exp(-value))

    def math_cosh(value):
        return 0.5 * (math.exp(value) + math.exp(-value))

    axis = 6378137.0  # GRS 80.
    flattening = 1.0 / 298.257222101  # GRS 80.

    central_meridian = 15.00
    scale = 0.9996
    false_northing = 0.0
    false_easting = 500000.0

    lon_lat = [0] * 2

    # Prepare ellipsoid-based stuff.
    e2 = flattening * (2.0 - flattening)
    n = flattening / (2.0 - flattening)
    a_roof = axis / (1.0 + n) * (1.0 + n*n/4.0 + n*n*n*n/64.0)
    delta1 = n/2.0 - 2.0*n*n/3.0 + 37.0*n*n*n/96.0 - n*n*n*n/360.0
    delta2 = n*n/48.0 + n*n*n/15.0 - 437.0*n*n*n*n/1440.0
    delta3 = 17.0*n*n*n/480.0 - 37*n*n*n*n/840.0
    delta4 = 4397.0*n*n*n*n/161280.0

    Astar = e2 + e2*e2 + e2*e2*e2 + e2*e2*e2*e2
    Bstar = -(7.0*e2*e2 + 17.0*e2*e2*e2 + 30.0*e2*e2*e2*e2) / 6.0
    Cstar = (224.0*e2*e2*e2 + 889.0*e2*e2*e2*e2) / 120.0
    Dstar = -(4279.0*e2*e2*e2*e2) / 1260.0

    # Convert.
    deg_to_rad = math.pi / 180
    lambda_zero = central_meridian * deg_to_rad
    xi = (x - false_northing) / (scale * a_roof)
    eta = (y - false_easting) / (scale * a_roof)
    xi_prim = xi - delta1*math.sin(2.0*xi) * math_cosh(2.0*eta) - delta2*math.sin(4.0*xi) * math_cosh(4.0*eta) - delta3*math.sin(6.0*xi) * math_cosh(6.0*eta) - delta4*math.sin(8.0*xi) * math_cosh(8.0*eta)
    eta_prim = eta - delta1*math.cos(2.0*xi) * math_sinh(2.0*eta) - delta2*math.cos(4.0*xi) * math_sinh(4.0*eta) - delta3*math.cos(6.0*xi) * math_sinh(6.0*eta) - delta4*math.cos(8.0*xi) * math_sinh(8.0*eta)
    phi_star = math.asin(math.sin(xi_prim) / math_cosh(eta_prim))
    delta_lambda = math.atan(math_sinh(eta_prim) / math.cos(xi_prim))
    lon_radian = lambda_zero + delta_lambda
    lat_radian = phi_star + math.sin(phi_star) * math.cos(phi_star) * (Astar + Bstar*math.pow(math.sin(phi_star), 2) + Cstar*math.pow(math.sin(phi_star), 4) + Dstar*math.pow(math.sin(phi_star), 6))
    lon_lat[1] = lat_radian * 180.0 / math.pi
    lon_lat[0] = lon_radian * 180.0 / math.pi
    return lon_lat

@st.cache_resource
def load_gp_data(filename, layername):
    """
    Read geodata from geopckage file of polygons with coordinates
    expressed as SWEREF 99 TM, with x=E y=N
    and transform to lat, lon.

    @param filename: name and location of geopackage file
    @param layername: name of layer to load
    @return: a dataframe with the coords expressed as lat, lon
    """
    gdata = gpd.read_file(filename, layers=layername)

    transform_df_coords(gdata, grid_to_geodetic)

    df_new = pd.DataFrame()
    df_new["coordinates"] = gdata["geometry"].apply(lambda g: [[list(t) for t in g.__geo_interface__['coordinates'][0]]])
    # df_new["elevation"] = np.log(gdata["Pop"]+1)*1000
    df_new["elevation"] = gdata["Pop"]+100

    return df_new


def gmap_widget(data, lat, lon, zoom):
    st.write(
        pdk.Deck(
            map_style="mapbox://styles/mapbox/light-v9",
            initial_view_state={
                "latitude": lat,
                "longitude": lon,
                "zoom": zoom,
                "pitch": 60,
            },
            layers=[
                pdk.Layer(
                    "PolygonLayer",
                    data,
                    id="geojson",
                    opacity=0.5,
                    stroked=False,
                    get_polygon="coordinates",
                    filled=True,
                    extruded=True,
                    wireframe=True,
                    get_elevation="elevation",
                    get_fill_color=[255, 255, 100],
                    get_line_color=[255, 255, 255],
                    auto_highlight=True,
                    pickable=False,
                )
                ,
            ],

        )
    )

gdata = load_gp_data(filename="~/Data/SCB/totalbefolkning_1km_211231/Totalbefolkning_1km_211231.gpkg", layername="Totalbefolkning_1km_211231")

row1_1, row1_2 = st.columns((2, 3))

with row1_1:
    st.title("Statistics Sweden open geodata - population")

with row1_2:
    st.write(
        """
    ##
    Test
    """
    )

map_container = st.container()

uppsala = [59.86109193413548, 17.627849966119896]
sydsverige = [56.249417938658254, 15.597110073971624]

with map_container:
    gmap_widget(gdata, uppsala[0], uppsala[1], 9)
