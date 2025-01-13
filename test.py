from numpy import inf
import networkx as nx
import osmnx as ox
import copy
import matplotlib.pyplot as plt
from datetime import datetime
import itertools
import statistics

def lex_max_flow(netzwerk, reihenfolge, zeitgrenze):
    #supply/demand-Vektor initialisieren: 
    x = [0 for i in reihenfolge]
    #counter anlegen
    j = 0
    #netzwerk zum expanded netzwerk formen:
    for i in reihenfolge:
        if netzwerk.nodes[i]['terminal'] == 1:
            netzwerk.add_edge(0, i, 0, capacity = inf, weight = 0)
    #jetzt der produktive Teil vom Algorithmus
    for i in reihenfolge:
        if netzwerk.nodes[i]['terminal'] == -1:
            #kante in netzwerk hinzufügen: 
            netzwerk.add_edge(i, 0, 0, capacity = inf, weight = -zeitgrenze)
            #fancy built_in min-cost-flow-finder
            y_costs, y_flow  = nx.network_simplex(netzwerk)
            #manipuliere netzwerk entsprechend dem Flow
            residual_graph_anlegen(y_flow, netzwerk)
        if netzwerk.nodes[i]['terminal'] == 1:
            #vorbereiten, dass die Kante von 0 zu i flow null haben soll:
            #zu diesem zwecke: gib 0 einen Lieferwert der dem kumulierten flow auf der Kante (0, i) entspricht und i einen entsprechenden Bedarf
            #der kumulierte Flow ist kapazität der residual-kante
            #i zu 0 ist eine backward arc und hat daher den key -1
            #teste, ob überhaupt was zu tun ist
            try:
                #supply von 0 setzen
                netzwerk.nodes[0]['demand'] = -netzwerk[i][0][-1]['capacity']
                #demand von i setzen
                netzwerk.nodes[i]['demand'] = netzwerk[i][0][-1]['capacity']
                #hier also keine min-cost-circ sondern min-cost-flow mit supply und Demand, wo die Kante (0, i) nicht genutzt werden darf
                netzwerk.remove_edge(0, i, 0)
                y_costs, y_flow  = nx.network_simplex(netzwerk)
                #manipuliere netzwerk entsprechend dem Flow
                residual_graph_anlegen(y_flow, netzwerk)
                #y_cost ist nicht der korrekte Wert, denn es fehlt die Kante von i zu 0, welche die komplette Kapazität ausschöpft
                y_costs += netzwerk[i][0][-1]['capacity'] * netzwerk[i][0][-1]['weight']
                #y_flow ist nicht der vollständige flow, es fehlt noch der Flow auf der Kante von i zu 0
                y_flow[i].update({0: {-1: netzwerk.nodes[i]['demand']}})
                #kante abziehen:
                netzwerk.remove_edge(i, 0, -1)
                #supply und demand zurücksetzen
                netzwerk.nodes[0]['demand'] = 0
                netzwerk.nodes[i]['demand'] = 0
            except KeyError:
                netzwerk.remove_edge(0, i, 0)
                y_costs = 0
        #supply/demand-Vektor aktualisieren
        #für den i-ten Reihenfolge-Knoten sind das: (kosten von flow yi) - (kosten von yi-1)
        x[j] = y_costs
        j += 1
        #füge der ausgabe den supply/demand-vektor hinzu
    print(sum(x))
    return x#,gamma

def residual_graph_anlegen(flow, G):
    for u in flow:
        for v in flow[u]:
            for i in flow[u][v]:
                if flow[u][v][i] > 0:
                    #1. capacity der arc reduzieren, die genutzt wurde
                    G[u][v][i]['capacity'] -= flow[u][v][i]
                    #2. capacity der backward-arc erhöhen (und anlegen falls nicht existent)
                    try:
                        G[v][u][-i-1]['capacity'] += flow[u][v][i]
                    except KeyError:
                        #backward arcs anlegen, falls noch nicht existent:
                        G.add_edge(v, u, -i-1, capacity = flow[u][v][i], weight = -G[u][v][i]['weight'], flow = 0)

#!!!!DAS DING FUNKTIONIERT FÜR BELIEBIGE KANTENZAHLEN!!!!
def flow_decomp(flow, G):
    #flow als Kanteneigenschaft setzen
    for u in flow:
        for v in flow[u]:
            for i in flow[u][v]:
                G[u][v][i]['flow'] = flow[u][v][i]
    #netzwerk duplizieren, um es leichter zu manipulieren
    netzwerk = copy.deepcopy(G)
    pfade = []
    #solange es noch kanten aus 0 gibt mit positiven Flow beginne loszulaufen
    while len(list(netzwerk.successors(0))) > 0:
        #laufe in 0 los
        pfad = [(0, 0)]
        min_flow = inf
        #while nicht wieder in 0 gehe positiv-flow Kanten entlang
        while pfad[-1][0] != 0 or pfad == [(0, 0)]:
            #wähle den nächsten knoten als "beliebigen" Nachfolger
            next_node = list(netzwerk.successors(pfad[-1][0]))[0]
            for key in netzwerk[pfad[-1][0]][next_node]:
                if netzwerk[pfad[-1][0]][next_node][key]['flow'] == 0:
                    #wenn kein flow auf der kante ist lösche die kante
                    netzwerk.remove_edge(pfad[-1][0], next_node, key)
                    #wenn es jetzt keine Kante mehr aus 0 rausgibt, dann sind wir hier fertig
                    if list(netzwerk.successors(pfad[-1][0])) == []:
                        return pfade
                    break
                #wenn sich ein teilkreis schließt, dann muss der flow auf dem kreis reduziert werden, dann kann weitergelaufen werden, bevor der kreis begonnen hat
                elif next_node in [k[0] for k in pfad] and next_node != 0:
                    min_flow_sub_circ = inf
                    #finde den bottleneck des teilkreises
                    for i in range(pfad.index((next_node, key)), len(pfad)-1):
                        min_flow_sub_circ = min(netzwerk[pfad[i]][pfad[i+1]][key]['flow'], min_flow_sub_circ)
                    #reduziere den flow auf dem teilkreis und lösche die bottleneck arc
                    for i in range(pfad.index((next_node, key)), len(pfad)-1):
                        netzwerk[pfad[i]][pfad[i+1]][key]['flow'] -= min_flow_sub_circ
                        if netzwerk[pfad[i]][pfad[i+1]][key]['flow'] == 0:
                           netzwerk.remove_edge(pfad[i], pfad[i+1], key) 
                    #setze den pfad zurück
                    pfad = pfad[:pfad.index((next_node, key))+1]
                    break
                else:
                    #wenn alles läuft wie gewünscht
                    min_flow = min(min_flow, netzwerk[pfad[-1][0]][next_node][key]['flow'])
                    pfad.append((next_node, key))
                    break
        #reduziere den flow auf dem teilkreis und lösche die bottleneck arc
        for i in range(len(pfad)-1):
            netzwerk[pfad[i][0]][pfad[i+1][0]][pfad[i+1][1]]['flow'] -= min_flow
            if netzwerk[pfad[i][0]][pfad[i+1][0]][pfad[i+1][1]]['flow'] == 0:
                netzwerk.remove_edge(pfad[i][0], pfad[i+1][0], pfad[i+1][1])
        pfade.append([pfad, min_flow])
    return pfade

def capacities_und_weight_setzen_osm_graph(Graph):
    Graph = ox.routing.add_edge_speeds(Graph)
    Graph = ox.routing.add_edge_travel_times(Graph)
    for u in Graph.nodes:
        for v in Graph[u]:
            for i in Graph[u][v]:
                Graph[u][v][i]['weight'] = int(Graph[u][v][i]['travel_time']*10)
                try:
                    Graph[u][v][i]['capacity'] = int(min(Graph[u][v][i]['lanes']))
                except:
                    Graph[u][v][i]['capacity'] = 1
                try:
                    del Graph[u][v][i]['speed_kph']
                    del Graph[u][v][i]['name']
                    del Graph[u][v][i]['highway']
                    del Graph[u][v][i]['maxspeed']
                    del Graph[u][v][i]['oneway']
                    del Graph[u][v][i]['reversed']
                    del Graph[u][v][i]['length']
                    del Graph[u][v][i]['travel_time']
                    del Graph[u][v][i]['geometry']
                    del Graph[u][v][i]['lanes']
                except:
                    pass


    
def mache_sekunden_mal_10_aus_meiner_minutenangabe(minuten):
    return minuten*600

def alle_permutationen(liste):
    return [list(p) for p in itertools.permutations(liste)]

#1 liefert
#-1 fordert


#print("Startzeit vom Graph=", datetime.now().strftime("%H:%M:%S"))
#G = ox.io.load_graphml('polen.graphml')
#print("Endzeit vom Graphen=", datetime.now().strftime("%H:%M:%S"))

#Graphen laden und speichern
'''
laender_liste=["Phoenix", "Philadelphia", "San Diego"]

for ort in laender_liste:
    try:
        print("Start " + ort)
        print("Startzeit " + ort + "importieren=", datetime.now().strftime("%H:%M:%S.%f"))
        G = ox.graph_from_place(ort, network_type="drive", simplify=True)
        #print("Endzeit " + ort + " importieren=", datetime.now().strftime("%H:%M:%S"))
        
        print("Anzahl Kanten " + ort + ":" + str(len(list(G.edges))))
        print("Anzahl Knoten " + ort + ":" + str(len(list(G.nodes))))
        #ort    kanten  knoten
        print(ort + "   " + str(len(list(G.edges)))+ "   " + str(len(list(G.nodes))))
        #print("Startzeit speichern=", datetime.now().strftime("%H:%M:%S"))
        ox.io.save_graphml(G, filepath=ort+".graphml", gephi=False, encoding='utf-8')
        print("Endzeit speichern=", datetime.now().strftime("%H:%M:%S"))
        print("Ende " + ort)
    except:
        print(ort + " wurde nicht gefunden")
'''


def terminals_setzen_reihenfolge_ausgeben(G, koordinaten_liste):
    reihenfolge=[]
    for koordinate in koordinaten_liste:
        zur_koordinate_passender_knoten = ox.nearest_nodes(G, koordinate[1], koordinate[0])
        G.nodes[zur_koordinate_passender_knoten]['terminal'] = koordinate[2]
        if koordinate[2]==1:
            reihenfolge.append(zur_koordinate_passender_knoten)
        else:
            reihenfolge.insert(0, zur_koordinate_passender_knoten)
    return reihenfolge


def laufzeit_auswertung(ortsauswahl):
    laufzeit_liste=[[], [], []]
    #wieso kann der Graph nicht hier gespeichert und immer wieder aufgerufen werden??????????
    for ort in ortsauswahl:
        print("Startzeit " + ort, datetime.now().strftime("%H:%M:%S.%f"))
        for i in [5, 10, 15]:
            G = ox.io.load_graphml(ort + ".graphml")
            reihenfolge = terminals_setzen_reihenfolge_ausgeben(G, orte[ort]["koordinaten_liste"][:i])
            capacities_und_weight_setzen_osm_graph(G)
            zeitgrenze = mache_sekunden_mal_10_aus_meiner_minutenangabe(600)
            print("Startzeit " + ort + " mit " + str(i) + " Terminalen", datetime.now().strftime("%H:%M:%S.%f"))
            startzeit=datetime.now()
            print(lex_max_flow(G, reihenfolge, zeitgrenze))
            endzeit = datetime.now()
            print("Endzeit " + ort + " mit " + str(i) + " Terminalen", datetime.now().strftime("%H:%M:%S.%f"))
            dauer = endzeit - startzeit
            print("Dauer: ", dauer)
            print(dauer.seconds)
            print(dauer.microseconds)
            dauer_in_microseconds = dauer.seconds * 1000000 + dauer.microseconds
            print(dauer_in_microseconds)
            laufzeit_liste[int((i/5)-1)].append(dauer_in_microseconds)
        print("Endzeit " + ort, datetime.now().strftime("%H:%M:%S.%f"))
    laufzeit_5_terminale_10k = laufzeit_liste[0][:3]
    laufzeit_5_terminale_100k = laufzeit_liste[0][3:6]
    laufzeit_5_terminale_1m = laufzeit_liste[0][6:]
    laufzeit_10_terminale_10k = laufzeit_liste[1][:3]
    laufzeit_10_terminale_100k = laufzeit_liste[1][3:6]
    laufzeit_10_terminale_1m = laufzeit_liste[1][6:]
    laufzeit_15_terminale_10k = laufzeit_liste[2][:3]
    laufzeit_15_terminale_100k = laufzeit_liste[2][3:6]
    laufzeit_15_terminale_1m = laufzeit_liste[2][6:]
    durchschnitt_laufzeit_5_terminale_10k = int(statistics.mean(laufzeit_5_terminale_10k))/1000000
    durchschnitt_laufzeit_5_terminale_100k = int(statistics.mean(laufzeit_5_terminale_100k))/1000000
    durchschnitt_laufzeit_5_terminale_1m = int(statistics.mean(laufzeit_5_terminale_1m))/1000000
    durchschnitt_laufzeit_10_terminale_10k = int(statistics.mean(laufzeit_10_terminale_10k))/1000000
    durchschnitt_laufzeit_10_terminale_100k = int(statistics.mean(laufzeit_10_terminale_100k))/1000000
    durchschnitt_laufzeit_10_terminale_1m = int(statistics.mean(laufzeit_10_terminale_1m))/1000000
    durchschnitt_laufzeit_15_terminale_10k = int(statistics.mean(laufzeit_15_terminale_10k))/1000000
    durchschnitt_laufzeit_15_terminale_100k = int(statistics.mean(laufzeit_15_terminale_100k))/1000000
    durchschnitt_laufzeit_15_terminale_1m = int(statistics.mean(laufzeit_15_terminale_1m))/1000000
    intervall_laufzeit_5_terminale_10k = (min(laufzeit_5_terminale_10k), durchschnitt_laufzeit_5_terminale_10k, max(laufzeit_5_terminale_10k), laufzeit_5_terminale_10k)
    intervall_laufzeit_5_terminale_100k = (min(laufzeit_5_terminale_100k), durchschnitt_laufzeit_5_terminale_100k, max(laufzeit_5_terminale_100k), laufzeit_5_terminale_100k)
    intervall_laufzeit_5_terminale_1m = (min(laufzeit_5_terminale_1m), durchschnitt_laufzeit_5_terminale_1m, max(laufzeit_5_terminale_1m), laufzeit_5_terminale_1m)
    intervall_laufzeit_10_terminale_10k = (min(laufzeit_10_terminale_10k), durchschnitt_laufzeit_10_terminale_10k, max(laufzeit_10_terminale_10k), laufzeit_10_terminale_10k)
    intervall_laufzeit_10_terminale_100k = (min(laufzeit_10_terminale_100k), durchschnitt_laufzeit_10_terminale_100k, max(laufzeit_10_terminale_100k), laufzeit_10_terminale_100k)
    intervall_laufzeit_10_terminale_1m = (min(laufzeit_10_terminale_1m), durchschnitt_laufzeit_10_terminale_1m, max(laufzeit_10_terminale_1m), laufzeit_10_terminale_1m)
    intervall_laufzeit_15_terminale_10k = (min(laufzeit_15_terminale_10k), durchschnitt_laufzeit_15_terminale_10k, max(laufzeit_15_terminale_10k), laufzeit_15_terminale_10k)
    intervall_laufzeit_15_terminale_100k = (min(laufzeit_15_terminale_100k), durchschnitt_laufzeit_15_terminale_100k, max(laufzeit_15_terminale_100k), laufzeit_15_terminale_100k)
    intervall_laufzeit_15_terminale_1m = (min(laufzeit_15_terminale_1m), durchschnitt_laufzeit_15_terminale_1m, max(laufzeit_15_terminale_1m), laufzeit_15_terminale_1m)
    return "intervall_laufzeit_5_terminale_10k: " + str(intervall_laufzeit_5_terminale_10k) + "\n" + "intervall_laufzeit_5_terminale_100k: " + str(intervall_laufzeit_5_terminale_100k) + "\n" + "intervall_laufzeit_5_terminale_1m: " + str(intervall_laufzeit_5_terminale_1m) + "\n" + "intervall_laufzeit_10_terminale_10k: " + str(intervall_laufzeit_10_terminale_10k) + "\n" + "intervall_laufzeit_10_terminale_100k: " + str(intervall_laufzeit_10_terminale_100k) + "\n" + "intervall_laufzeit_10_terminale_1m: " + str(intervall_laufzeit_10_terminale_1m) + "\n" + "intervall_laufzeit_15_terminale_10k: " + str(intervall_laufzeit_15_terminale_10k) + "\n" + "intervall_laufzeit_15_terminale_100k: " + str(intervall_laufzeit_15_terminale_100k) + "\n" + "intervall_laufzeit_15_terminale_1m: " + str(intervall_laufzeit_15_terminale_1m)


orte = {
    "Brügge": {
        "koordinaten_liste": [
            (51.223205, 3.313496, 1),
            (51.338056, 3.229773, -1),
            (51.339499, 3.179540, -1),
            (51.239475, 3.161063, -1),
            (51.160241, 3.132770, 1),
            (51.250680, 3.250560, 1),
            (51.285722, 3.275965, -1),
            (51.178705, 3.233815, 1),
            (51.188115, 3.288090, -1),
            (51.163862, 3.212451, 1),
            (51.287167, 3.174343, 1),
            (51.241644, 3.266149, -1),
            (51.270914, 3.263262, -1),
            (51.206568, 3.147205, 1),
            (51.298360, 3.260664, 1)
        ]
    },
    "Manhattan": {
        "koordinaten_liste": [
            (40.878734, -73.921875, 1),
            (40.799849, -73.918099, 1),
            (40.747528, -73.966413, 1),
            (40.700707, -74.015203, -1),
            (40.726966, -74.016069, -1),
            (40.845575, -73.924622, 1),
            (40.714275, -73.973341, -1),
            (40.707491, -73.997881, -1),
            (40.761525, -74.007119, -1),
            (40.850464, -73.947070, 1),
            (40.877176, -73.905739, 1),
            (40.755265, -73.949425, -1),
            (40.828200, -73.932826, -1),
            (40.831397, -73.951230, -1),
            (40.735551, -73.975042, 1)
        ]
    },
    "Kassel": {
        "koordinaten_liste": [
            (51.346313, 9.475414, 1),
            (51.342453, 9.519359, 1),
            (51.275595, 9.496507, -1),
            (51.269725, 9.489990, -1),
            (51.265480, 9.451954, -1),
            (51.352723, 9.555612, 1),
            (51.267512, 9.434416, -1),
            (51.279252, 9.404968, -1),
            (51.275898, 9.529017, 1),
            (51.312951, 9.355485, -1),
            (51.314651, 9.539280, 1),
            (51.311908, 9.560409, 1),
            (51.293401, 9.559796, 1),
            (51.342228, 9.372357, -1),
            (51.349849, 9.416202, -1)
        ]
    },
    "Indiana": {
        "koordinaten_liste": [
            (41.779186, -86.729170, 1),
            (37.746842, -87.085476, 1),
            (41.155526, -84.759902, 1),
            (39.431244, -87.614312, -1),
            (40.781341, -87.587379, -1),
            (41.798668, -84.978458, 1),
            (41.708848, -87.531629, -1),
            (38.230330, -88.058523, -1),
            (38.230795, -85.687518, 1),
            (37.829443, -87.595011, -1),
            (40.121087, -87.572998, -1),
            (41.570841, -87.594652, -1),
            (41.642641, -84.743287, 1),
            (39.844211, -84.772967, 1),
            (39.249824, -84.772967, 1)
        ]
    },
    "Tennessee": {
        "koordinaten_liste": [
            (35.919785, -82.541698, 1),
            (35.165008, -90.197799, -1),
            (34.961824, -89.988759, 1),
            (36.635982, -82.202007, -1),
            (36.719808, -86.526529, -1),
            (35.729106, -83.038169, 1),
            (36.183857, -89.701328, -1),
            (36.531072, -88.891297, -1),
            (36.635982, -82.580893, -1),
            (34.940406, -86.905414, 1),
            (36.635982, -84.109500, -1),
            (36.677907, -87.375755, -1),
            (34.951116, -85.154702, 1),
            (34.940406, -85.468262, 1),
            (34.972530, -89.649068, 1)
        ]
    },
    "Tschechien": {
        "koordinaten_liste": [
            (50.840150, 15.512241, 1),
            (48.628245, 17.014718, 1),
            (49.965803, 18.478000, 1),
            (49.636930, 12.468092, -1),
            (48.800654, 15.041900, -1),
            (48.938156, 17.929270, 1),
            (50.425806, 16.296142, 1),
            (50.840150, 13.892179, -1),
            (48.576407, 14.467040, -1),
            (49.788995, 18.687041, 1),
            (50.707958, 16.008712, 1),
            (49.373924, 18.464935, 1),
            (48.714524, 16.635833, -1),
            (48.731762, 16.060972, -1),
            (49.092394, 13.186668, -1)
        ]
    },
    "San Antonio": {
        "koordinaten_liste": [
            (29.710963, -98.449899, 1),
            (29.488378, -98.222894, 1),
            (29.174832, -98.429727, 1),
            (29.252442, -98.585624, -1),
            (29.558397, -98.673231, -1),
            (29.577542, -98.338846, 1),
            (29.300792, -98.369677, 1),
            (29.289713, -98.642209, -1),
            (29.484144, -98.718346, -1),
            (29.656818, -98.521220, -1),
            (29.554814, -98.343745, 1),
            (29.391665, -98.353544, 1),
            (29.247404, -98.473609, 1),
            (29.389374, -98.673388, -1),
            (29.684266, -98.639341, -1)
        ]
    },
    "Fortaleza": {
        "koordinaten_liste": [
            (-3.721651, -38.509445, 1),
            (-3.811580, -38.429447, 1),
            (-3.888571, -38.501713, 1),
            (-3.840911, -38.584390, -1),
            (-3.761240, -38.613636, -1),
            (-3.707138, -38.470759, 1),
            (-3.873296, -38.489668, 1),
            (-3.821153, -38.618481, -1),
            (-3.748049, -38.601163, -1),
            (-3.699478, -38.588616, -1),
            (-3.814839, -38.423731, 1),
            (-3.854965, -38.472929, 1),
            (-3.859242, -38.542133, 1),
            (-3.776569, -38.621891, -1),
            (-3.740862, -38.601781, -1)
        ]
    },
    "Casablanca": {
        "koordinaten_liste": [
            (33.646792, -7.485928, 1),
            (33.625181, -7.465467, 1),
            (33.563870, -7.496831, 1),
            (33.527627, -7.565052, 1),
            (33.517494, -7.695443, -1),
            (33.640442, -7.480640, 1),
            (33.583136, -7.477288, 1),
            (33.509308, -7.611915, -1),
            (33.503123, -7.681970, -1),
            (33.565652, -7.712922, -1),
            (33.631346, -7.472628, 1),
            (33.599403, -7.469533, 1),
            (33.539339, -7.540528, 1),
            (33.501006, -7.628395, -1),
            (33.546230, -7.704790, -1)
        ]
    }
}





geeignete_graphenliste = ["Brügge", "Kassel", "Manhattan", "San Antonio", "Fortaleza", "Casablanca", "Indiana", "Tennessee", "Tschechien"]


print(laufzeit_auswertung(geeignete_graphenliste))


'''
for ort in geeignete_graphenliste:
    print("Startzeit " + ort, datetime.now().strftime("%H:%M:%S.%f"))
    for i in [5, 10, 15]:
        G = ox.io.load_graphml(ort + ".graphml")
        reihenfolge = terminals_setzen_reihenfolge_ausgeben(G, orte[ort]["koordinaten_liste"][:i])
        capacities_und_weight_setzen_osm_graph(G)
        zeitgrenze = mache_sekunden_mal_10_aus_meiner_minutenangabe(600)
        print("Startzeit " + ort + " mit " + str(i) + " Terminalen", datetime.now().strftime("%H:%M:%S.%f"))
        startzeit=datetime.now()
        print(lex_max_flow(G, reihenfolge, zeitgrenze))
        endzeit = datetime.now()
        print("Endzeit " + ort + " mit " + str(i) + " Terminalen", datetime.now().strftime("%H:%M:%S.%f"))
        print("Dauer: ", startzeit-endzeit)
    print("Endzeit " + ort, datetime.now().strftime("%H:%M:%S.%f"))

'''



#capacities_und_weight_setzen_osm_graph(G)
#print(len(list(G.edges)))
#print(len(list(G.nodes)))

#54.774934,9.423078  flensburgx
#orig_node = ox.nearest_nodes(G, 54.774934, 9.423078)

#dest_node = ox.nearest_nodes(G, lat, lon)

#23.0429451, 54.3409021 #polen oben rechts
#22.8532551, 49.0246205 #polen unten rechts
#14.8974501, 50.9004874 #polen unten links
#14.4703971, 53.4204332 #polen oben links

#nrw: 
#6.3863601,50.4014374 unten links
#8.1153281,50.7125854 unten rechts
#8.5726391,52.0084803 oben rechts 
#6.9185391,52.1719543 oben links

#G.nodes[ox.nearest_nodes(G, 23.0429451, 54.3409021)]['terminal'] = -1
#G.nodes[ox.nearest_nodes(G, 22.8532551, 49.0246205)]['terminal'] = -1
#G.nodes[ox.nearest_nodes(G, 14.8974501, 50.9004874)]['terminal'] = 1
#G.nodes[ox.nearest_nodes(G, 14.4703971, 53.4204332)]['terminal'] = 1


#print("Startzeit =", datetime.now().strftime("%H:%M:%S"))
#print(lex_max_flow(G, [ox.nearest_nodes(G, 23.0429451, 54.3409021), ox.nearest_nodes(G, 22.8532551, 49.0246205), ox.nearest_nodes(G, 14.8974501, 50.9004874), ox.nearest_nodes(G, 14.4703971, 53.4204332)], mache_sekunden_mal_10_aus_meiner_minutenangabe(420)))
#print("Endzeit =", datetime.now().strftime("%H:%M:%S"))


'''
G.nodes[1419616355]['terminal'] = -1 #unten rechts 
G.nodes[51884259]['terminal'] = 1 #oben rechts
G.nodes[1688192973]['terminal'] = -1 #unten links
G.nodes[36552365]['terminal'] = 1 #oben links

print(lex_max_flow(G, [1419616355, 1688192973, 36552365, 51884259], 200000))

'''


'''
G = nx.MultiDiGraph()
G.add_edge(1,2,capacity=4,weight=1)
G.add_edge(1,5,capacity=1,weight=1)
G.add_edge(2,3,capacity=2,weight=5)
G.add_edge(4,5,capacity=5,weight=4)
G.add_edge(5,6,capacity=7,weight=3)
G.add_edge(5,7,capacity=4,weight=1)
G.add_edge(7,5,capacity=8,weight=4)
G.add_edge(7,8,capacity=4,weight=1)
G.add_edge(8,3,capacity=3,weight=2)
G.nodes[3]['terminal'] = -1
G.nodes[6]['terminal'] = -1
G.nodes[1]['terminal'] = 1
G.nodes[4]['terminal'] = 1
for i in alle_permutationen([3,6,1,4]):
    print(i)
    print(lex_max_flow(G, i, 10))
'''
'''
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull

A = np.array([[14, -17,-12], [14, -29, 0], [14, -29, 0], [14, -29, 0], [14, -29, 0], [14, -29, 0], [14, 0, -29], [14, 0, -29], [14, 0, -29], [14, 0, -29], [14, 0, -29], [14, 0, -29], [0, -15, 0], [0, -15, 0], [0, 0, -15], [0, 0, -15], [0, 0, 0], [14, -14, 0], [14, -14, 0], [14, 0, -14], [14, 0, -14]])
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
for cube, color in zip([A], ['b']):
    hull = ConvexHull(cube)
    # draw the polygons of the convex hull
    for s in hull.simplices:
        tri = Poly3DCollection([cube[s]])
        tri.set_color(color)
        tri.set_alpha(0.5)
        ax.add_collection3d(tri)
    # draw the vertices
    ax.scatter(cube[:, 0], cube[:, 1], cube[:, 2], marker='o', color='purple')
plt.show()
'''




'''
1, 4, 7, 9, 13
2, 6, 10, 12, 15 
3, 5, 8, 11, 14

'''