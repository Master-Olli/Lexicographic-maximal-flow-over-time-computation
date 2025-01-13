from numpy import inf
import networkx as nx
import osmnx as ox
import copy
import matplotlib.pyplot as plt
import itertools


def lex_max_flow(netzwerk, reihenfolge, zeitgrenze):
    '''
    Die Funktion erwartet:
    ein gerichtetes Netzwerk als networkx graph
    eine Liste von Knoten (in Form von ihrer Knoten-id) in Reihenfolge, nach der sie lexikographisch maximiert werden
    eine Zeitgrenze als ganzzahlige Zahl

    Output ist:
    ein Vektor x, der den Verbrauch der Knoten angibt
    eine Pfadzerlegung gamma, nach deren Muster ein lexikographisch maximaler flow over time definiert werden kann.
    '''
    #supply/demand-Vektor initialisieren: 
    x = [0 for i in reihenfolge]
    #path-decomp initialisieren
    gamma = []
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
            #flow dekomponieren und chain-decomposition erweitern
            gamma.append(flow_decomp(y_flow, netzwerk))
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
                #flow dekomponieren und chain-decomposition erweitern
                gamma.append(flow_decomp(y_flow, netzwerk))
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
    return x, gamma

def residual_graph_anlegen(flow, G):
    #für alle kanten:
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
                #travel_time ist immer gesetzt, muss aber ganzzahlig werden, weswegen die Eingabe-Zeit mit 10 multipliziert werden muss
                Graph[u][v][i]['weight'] = int(Graph[u][v][i]['travel_time']*10)
                #falls die lanes als Eigenschaft nicht gesetzt sind, dann setze sie default auf 1
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
    #notwendig, da die travel_time unter Umständen nicht ganzzahlig ist sonder eine Nachkommastelle hat
    return minuten*600



#Beispielaufruf der Berechnung auf der Straßenkarte von "Charlottenburg-Wilmersdorf" in Berlin
G = ox.graph_from_place("Charlottenburg-Wilmersdorf", network_type="drive", simplify=True)
capacities_und_weight_setzen_osm_graph(G)

#Terminale setzen für den Beispielaufruf
G.nodes[484245]['terminal'] = 1
G.nodes[251105127]['terminal'] = 1
G.nodes[28200228]['terminal'] = -1
G.nodes[936149309]['terminal'] = 1
G.nodes[247078355]['terminal'] = -1
print(lex_max_flow(G, [247078355, 936149309, 28200228, 251105127, 484245], 9000))




#weitere Beispielaufrufe

G = ox.graph_from_place("Berlin", network_type="drive", simplify=True)
capacities_und_weight_setzen_osm_graph(G)

G.nodes[386026173]['terminal'] = 1 #unten rechts Köpenick
G.nodes[249812374]['terminal'] = -1 #oben links Spandau
G.nodes[1602528344]['terminal'] = -1 #unten links Wannsee
G.nodes[5412030559]['terminal'] = 1 #oben rechts weißensee
G.nodes[410940190]['terminal'] = -1 #unten mitte Mariendorf


G = ox.graph_from_place("Brandenburg", network_type="drive", simplify=True)
capacities_und_weight_setzen_osm_graph(G)

G.nodes[1419616355]['terminal'] = -1 #unten rechts 
G.nodes[51884259]['terminal'] = 1 #oben rechts
G.nodes[1688192973]['terminal'] = -1 #unten links
G.nodes[36552365]['terminal'] = 1 #oben links
print(len(list(G.edges)))
print(len(list(G.nodes)))

print(lex_max_flow(G, [1419616355, 1688192973, 36552365, 51884259], 200000))
