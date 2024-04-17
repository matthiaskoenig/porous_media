# load mesh with geojson
# pip install geojson descartes

import geojson
from pathlib import Path
import matplotlib.pyplot as plt
from descartes import PolygonPatch
from porous_media.console import console
import pandas as pd




path = Path("/home/mkoenig/work/qualiperf/zonation_meshes/image_data/control/LobuliSegmentation/MNT-021/0")
geojson_path = path / "lobuli.geojson"

distances_path = Path("/home/mkoenig/work/qualiperf/zonation_meshes/image_data/control/PortalityMap/lobule_distances.csv")
df = pd.read_csv(distances_path)
console.print(df.head())

df = df[(df.roi==0) & (df.subject=="MNT-021") & (df.protein=="cyp1a2")].copy()
df = df[(df.width>150) & (df.width<200) & (df.width>150) & (df.height<200)].copy()
print(len(df))

# x: width
# y: height
# intensity: protein amount (relative expression)
# pv_dist: portality


# ---





with open(geojson_path, "r") as f:
    data = geojson.load(f)

import matplotlib
fig = plt.figure(figsize=(20, 20), dpi=300)
ax = fig.gca()
BLUE = '#6699cc'

for kg, lobulus_geometry in enumerate(data["features"]):
    console.rule(style="white")
    console.print(lobulus_geometry)
    poly = lobulus_geometry["geometry"]

    # visualize geojson polygon
    ax.add_patch(PolygonPatch(poly, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2 ))

    if kg > 10:
        break

for k, row in df.iterrows():
    x, y, intensity = row["width"], row["height"], row["intensity"]
    console.print(x, y, intensity)
    # relative protein amount
    cmap = matplotlib.cm.get_cmap('Spectral')
    color = cmap(intensity)
    ax.plot(x, y, color=color, marker="o")


ax.axis('scaled')
plt.show()

# position and protein information;
