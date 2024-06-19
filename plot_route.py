# plot.py

import networkx as nx
import matplotlib.pyplot as plt


def plot_solution(G, pos, cable_vars, cables, xy_position, node_labels):
    S = G.copy()
    good_edges = []
    clean_nodes = set()

    for v in cable_vars:
        if v.x > 0 and v.varName.startswith('y'):
            nodes = v.varName.replace('y[', '').replace(']', '').split(',')
            good_edges.append((int(nodes[0]), int(nodes[1])))
            clean_nodes.update(nodes)

    edge_delete = [e for e in G.edges() if e not in good_edges]
    S.remove_edges_from(edge_delete)
    S.remove_nodes_from(set(G.nodes()) - clean_nodes)

    activated_edges = {cable['type']: [] for cable in cables}
    for i, j in S.edges():
        for cable in range(len(cables)):
            if cable_vars[i, j, cable].x > 0:
                activated_edges[cables[cable]['type']].append((i, j))

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_position([0.12, 0.11, 0.8, 0.8])
    colors = {'black': 'black', 'red': 'darkred', 'blue': 'darkblue'}
    widths = {'black': 2, 'red': 4, 'blue': 6}

    for color in colors:
        nx.draw_networkx_edges(S, pos, edgelist=activated_edges[color], edge_color=colors[color], width=widths[color],
                               arrows=False, ax=ax)

    ax.scatter(xy_position[:, 0], xy_position[:, 1], s=150, c='y', edgecolors='k')
    ax.grid(True, which='both', ls='-', linewidth=1, alpha=0.5)
    ax.set_xlabel('Cross-stream direction [m]')
    ax.set_ylabel('Stream-wise direction [m]')
    ax.set_aspect('equal')
    plt.show()
