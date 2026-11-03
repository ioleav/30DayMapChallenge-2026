"""
Red ferroviaria Maglev transcontinental (Cero Aviones) - 30 Day Map Challenge 2026

Genera cuatro archivos para Kepler.gl:
  - maglev_trips.geojson : trayectorias [lon, lat, 0, timestamp] -> capa Trip (luces animadas)
  - maglev_red.geojson   : la red estática (raíles) -> capa Line / GeoJSON
  - maglev_nodos.geojson : las metrópolis -> capa Point
  - maglev_arcos.csv     : pares origen-destino -> capa Arc (3D)
"""
import csv
import itertools
import json
import random

import geopandas as gpd
import networkx as nx
import numpy as np
from pyproj import Geod
from shapely.geometry import LineString, Point

# ------------------------------------------------------------------ parámetros
BASE_TS = 1700000000      # segundos Unix; Kepler necesita timestamps reales
VELOCIDAD_KMH = 1000      # velocidad media simulada del Maglev
PASO_KM = 100             # distancia entre vértices (curvatura geodésica y fluidez)
TRENES_POR_RUTA = 3       # salidas escalonadas por ruta
INTERVALO_H = 8           # horas entre salidas de una misma ruta
SEMILLA = 2026
geod = Geod(ellps="WGS84")

# ------------------------------------------------------------------ nodos (lon, lat, continente)
NODOS = {
    # Europa
    "Lisboa": (-9.14, 38.72, "Europa"), "Madrid": (-3.70, 40.42, "Europa"),
    "París": (2.35, 48.86, "Europa"), "Londres": (-0.13, 51.51, "Europa"),
    "Berlín": (13.40, 52.52, "Europa"), "Roma": (12.50, 41.90, "Europa"),
    "Estambul": (28.98, 41.01, "Europa"), "Moscú": (37.62, 55.75, "Europa"),
    # África
    "Tánger": (-5.80, 35.76, "África"), "Argel": (3.06, 36.75, "África"),
    "El Cairo": (31.24, 30.04, "África"), "Lagos": (3.38, 6.52, "África"),
    "Kinsasa": (15.27, -4.44, "África"), "Nairobi": (36.82, -1.29, "África"),
    "Yibuti": (43.15, 11.59, "África"), "Johannesburgo": (28.05, -26.20, "África"),
    "Ciudad del Cabo": (18.42, -33.92, "África"),
    # Asia
    "Saná": (44.19, 15.37, "Asia"), "Riad": (46.68, 24.71, "Asia"),
    "Teherán": (51.39, 35.69, "Asia"), "Dubái": (55.27, 25.20, "Asia"),
    "Delhi": (77.21, 28.61, "Asia"), "Bombay": (72.88, 19.08, "Asia"),
    "Daca": (90.41, 23.81, "Asia"), "Bangkok": (100.50, 13.76, "Asia"),
    "Singapur": (103.82, 1.35, "Asia"), "Yakarta": (106.85, -6.21, "Asia"),
    "Novosibirsk": (82.92, 55.03, "Asia"), "Ulán Bator": (106.92, 47.89, "Asia"),
    "Pekín": (116.41, 39.90, "Asia"), "Shanghái": (121.47, 31.23, "Asia"),
    "Seúl": (126.98, 37.57, "Asia"), "Tokio": (139.65, 35.68, "Asia"),
    "Yakutsk": (129.73, 62.03, "Asia"), "Uelen": (-169.82, 66.16, "Asia"),
    # América del Norte
    "Wales": (-168.09, 65.61, "América"), "Fairbanks": (-147.72, 64.84, "América"),
    "Vancouver": (-123.12, 49.28, "América"), "Los Ángeles": (-118.24, 34.05, "América"),
    "Chicago": (-87.63, 41.88, "América"), "Nueva York": (-74.01, 40.71, "América"),
    "Ciudad de México": (-99.13, 19.43, "América"),
    # América Central y del Sur
    "Panamá": (-79.52, 8.98, "América"), "Bogotá": (-74.07, 4.71, "América"),
    "Lima": (-77.04, -12.05, "América"), "São Paulo": (-46.63, -23.55, "América"),
    "Santiago": (-70.67, -33.45, "América"), "Buenos Aires": (-58.38, -34.60, "América"),
    # Oceanía
    "Darwin": (130.84, -12.46, "Oceanía"), "Brisbane": (153.03, -27.47, "Oceanía"),
    "Sídney": (151.21, -33.87, "Oceanía"), "Perth": (115.86, -31.95, "Oceanía"),
}

# ------------------------------------------------------------------ aristas (origen, destino, tipo)
# tipo "cruce" = túnel submarino en un punto crítico
ARISTAS = [
    # Europa
    ("Lisboa", "Madrid", "tierra"), ("Madrid", "París", "tierra"),
    ("París", "Londres", "cruce: Canal de la Mancha"), ("París", "Berlín", "tierra"),
    ("París", "Roma", "tierra"), ("Berlín", "Moscú", "tierra"),
    ("Roma", "Estambul", "tierra"), ("Berlín", "Estambul", "tierra"),
    ("Estambul", "Moscú", "tierra"),
    # Gibraltar y África
    ("Madrid", "Tánger", "cruce: Estrecho de Gibraltar"), ("Tánger", "Argel", "tierra"),
    ("Argel", "El Cairo", "tierra"), ("Argel", "Lagos", "tierra"),
    ("Lagos", "Kinsasa", "tierra"), ("Kinsasa", "Nairobi", "tierra"),
    ("Kinsasa", "Johannesburgo", "tierra"), ("Nairobi", "Johannesburgo", "tierra"),
    ("Johannesburgo", "Ciudad del Cabo", "tierra"), ("El Cairo", "Nairobi", "tierra"),
    ("El Cairo", "Yibuti", "tierra"), ("Nairobi", "Yibuti", "tierra"),
    ("Lagos", "El Cairo", "tierra"),
    # Oriente Medio y Suez
    ("El Cairo", "Estambul", "cruce: Canal de Suez / Anatolia"),
    ("Yibuti", "Saná", "cruce: Bab el-Mandeb"), ("Saná", "Riad", "tierra"),
    ("Riad", "Dubái", "tierra"), ("Riad", "Teherán", "tierra"),
    ("Teherán", "Estambul", "tierra"), ("Teherán", "Dubái", "cruce: Golfo Pérsico"),
    ("Dubái", "Bombay", "cruce: Mar Arábigo"), ("Teherán", "Delhi", "tierra"),
    # Asia
    ("Bombay", "Delhi", "tierra"), ("Delhi", "Daca", "tierra"),
    ("Daca", "Bangkok", "tierra"), ("Bangkok", "Singapur", "tierra"),
    ("Singapur", "Yakarta", "cruce: Estrecho de Sunda"),
    ("Delhi", "Novosibirsk", "tierra"), ("Moscú", "Novosibirsk", "tierra"),
    ("Novosibirsk", "Ulán Bator", "tierra"), ("Ulán Bator", "Pekín", "tierra"),
    ("Pekín", "Shanghái", "tierra"), ("Shanghái", "Bangkok", "tierra"),
    ("Pekín", "Seúl", "cruce: Mar Amarillo"), ("Seúl", "Tokio", "cruce: Estrecho de Corea"),
    ("Novosibirsk", "Yakutsk", "tierra"), ("Yakutsk", "Uelen", "tierra"),
    # Estrecho de Bering
    ("Uelen", "Wales", "cruce: Estrecho de Bering"), ("Wales", "Fairbanks", "tierra"),
    # América
    ("Fairbanks", "Vancouver", "tierra"), ("Vancouver", "Los Ángeles", "tierra"),
    ("Vancouver", "Chicago", "tierra"), ("Los Ángeles", "Chicago", "tierra"),
    ("Chicago", "Nueva York", "tierra"), ("Los Ángeles", "Ciudad de México", "tierra"),
    ("Chicago", "Ciudad de México", "tierra"), ("Ciudad de México", "Panamá", "tierra"),
    ("Panamá", "Bogotá", "tierra"), ("Bogotá", "Lima", "tierra"),
    ("Lima", "Santiago", "tierra"), ("Santiago", "Buenos Aires", "tierra"),
    ("Buenos Aires", "São Paulo", "tierra"), ("Bogotá", "São Paulo", "tierra"),
    ("Lima", "São Paulo", "tierra"),
    # Atlántico
    ("Lisboa", "Nueva York", "cruce: Atlántico Norte"),
    # Oceanía
    ("Yakarta", "Darwin", "cruce: Mar de Timor"), ("Darwin", "Brisbane", "tierra"),
    ("Brisbane", "Sídney", "tierra"), ("Darwin", "Perth", "tierra"),
    ("Perth", "Sídney", "tierra"),
]

# Metrópolis que actúan como orígenes/destinos de las rutas animadas
HUBS = [
    "Madrid", "Londres", "Moscú", "El Cairo", "Lagos", "Johannesburgo",
    "Dubái", "Delhi", "Singapur", "Pekín", "Tokio", "Sídney",
    "Vancouver", "Nueva York", "Ciudad de México", "São Paulo",
]


# ------------------------------------------------------------------ utilidades
def distancia_km(a, b):
    lon1, lat1 = NODOS[a][:2]
    lon2, lat2 = NODOS[b][:2]
    return geod.inv(lon1, lat1, lon2, lat2)[2] / 1000


def densificar(a, b):
    """Puntos [lon, lat] sobre el círculo máximo entre a y b (incluye extremos)."""
    lon1, lat1 = NODOS[a][:2]
    lon2, lat2 = NODOS[b][:2]
    n = max(1, int(distancia_km(a, b) // PASO_KM))
    medios = geod.npts(lon1, lat1, lon2, lat2, n) if n > 1 else []
    return [[lon1, lat1]] + [list(p) for p in medios] + [[lon2, lat2]]


def dividir_antimeridiano(puntos):
    """
    puntos: lista de [lon, lat, ...extra]. Parte la línea al cruzar lon=±180
    e interpola lat y el resto de valores (p. ej. el timestamp) en el corte.
    """
    tramos, actual = [], [puntos[0]]
    for prev, cur in zip(puntos[:-1], puntos[1:]):
        if abs(cur[0] - prev[0]) > 180:
            este = prev[0] > 0                       # cruza hacia +180 (este)
            lon_cur = cur[0] + 360 if este else cur[0] - 360
            borde = 180 if este else -180
            f = (borde - prev[0]) / (lon_cur - prev[0])
            extra = [p + f * (c - p) for p, c in zip(prev[1:], cur[1:])]
            actual.append([borde] + extra)
            tramos.append(actual)
            actual = [[-borde] + extra]
        actual.append(cur)
    tramos.append(actual)
    return tramos


# ------------------------------------------------------------------ grafo
G = nx.Graph()
for nombre, (lon, lat, cont) in NODOS.items():
    G.add_node(nombre, lon=lon, lat=lat, continente=cont)
for a, b, tipo in ARISTAS:
    G.add_edge(a, b, km=distancia_km(a, b), tipo=tipo)

assert nx.is_connected(G), "El grafo no está conectado"

# ------------------------------------------------------------------ GeoPandas: nodos y red estática
gdf_nodos = gpd.GeoDataFrame(
    [{"nombre": n, "continente": d["continente"], "grado": G.degree(n),
      "geometry": Point(d["lon"], d["lat"])} for n, d in G.nodes(data=True)],
    crs="EPSG:4326",
)

filas_red = []
for a, b, d in G.edges(data=True):
    for tramo in dividir_antimeridiano(densificar(a, b)):
        if len(tramo) > 1:
            filas_red.append({
                "origen": a, "destino": b, "km": round(d["km"]),
                "tipo": "cruce" if d["tipo"].startswith("cruce") else "tierra",
                "detalle": d["tipo"], "geometry": LineString(tramo),
            })
gdf_red = gpd.GeoDataFrame(filas_red, crs="EPSG:4326")

gdf_nodos.to_file("maglev_nodos.geojson", driver="GeoJSON")
gdf_red.to_file("maglev_red.geojson", driver="GeoJSON")

# ------------------------------------------------------------------ rutas animadas (trips)
random.seed(SEMILLA)
features, arcos = [], []
id_ruta = 0

for origen, destino in itertools.combinations(HUBS, 2):
    camino = nx.shortest_path(G, origen, destino, weight="km")

    puntos = []
    for a, b in zip(camino[:-1], camino[1:]):
        seg = densificar(a, b)
        puntos += seg if not puntos else seg[1:]
    lons, lats = np.array(puntos).T

    # Distancia acumulada -> tiempo a velocidad constante
    dist = geod.inv(lons[:-1], lats[:-1], lons[1:], lats[1:])[2] / 1000
    acum = np.concatenate([[0], np.cumsum(dist)])
    segundos = acum / VELOCIDAD_KMH * 3600
    km_total = float(acum[-1])
    horas = km_total / VELOCIDAD_KMH

    salida_base = random.uniform(0, INTERVALO_H * 3600)
    for tren in range(TRENES_POR_RUTA):
        t0 = BASE_TS + salida_base + tren * INTERVALO_H * 3600
        pts = [[float(lo), float(la), round(t0 + s)] for lo, la, s in zip(lons, lats, segundos)]
        for k, tramo in enumerate(dividir_antimeridiano(pts)):
            if len(tramo) < 2:
                continue
            features.append({
                "type": "Feature",
                "properties": {
                    "id_ruta": id_ruta, "tren": tren, "tramo": k,
                    "origen": origen, "destino": destino,
                    "continente_origen": NODOS[origen][2],
                    "km": round(km_total), "horas": round(horas, 1),
                },
                # [lon, lat, altitud, timestamp]
                "geometry": {"type": "LineString",
                             "coordinates": [[p[0], p[1], 0, p[2]] for p in tramo]},
            })
    arcos.append({
        "origen": origen, "destino": destino,
        "lon_origen": NODOS[origen][0], "lat_origen": NODOS[origen][1],
        "lon_destino": NODOS[destino][0], "lat_destino": NODOS[destino][1],
        "km": round(km_total), "continente_origen": NODOS[origen][2],
    })
    id_ruta += 1

with open("maglev_trips.geojson", "w", encoding="utf-8") as f:
    json.dump({"type": "FeatureCollection", "features": features}, f, ensure_ascii=False)

with open("maglev_arcos.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(arcos[0].keys()))
    w.writeheader()
    w.writerows(arcos)

print(f"Nodos: {G.number_of_nodes()} | Aristas: {G.number_of_edges()}")
print(f"Rutas: {id_ruta} | Features de trips: {len(features)}")
print(f"Duración total de la animación: ~{(max(f['geometry']['coordinates'][-1][3] for f in features) - BASE_TS) / 3600:.0f} h")
