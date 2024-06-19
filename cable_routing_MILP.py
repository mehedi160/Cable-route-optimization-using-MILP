import networkx as nx
import gurobipy as gb
import numpy as np
import math


class CableRoutingMILP:
    def __init__(self, layout_file, rows=10, cols=10, cell_width=50, cell_height=50, max_cable_subst=20):
        self.G = nx.Graph()
        self.cables = [
            {'type': "black", 'unit': 120, 'flow_max': 42},
            {'type': "red", 'unit': 150.5, 'flow_max': 48},
            {'type': "blue", 'unit': 182.4, 'flow_max': 53}
        ]
        self.max_cable_subst = max_cable_subst
        self.layout_file = layout_file
        self.rows = rows
        self.cols = cols
        self.cell_width = cell_width
        self.cell_height = cell_height
        self.pos = {}
        self.node_labels = {}
        self.xy_position, self.cr_position = self.initialize_positions()
        self.initialize_graph()
        self.create_optimization_model()

    def initialize_positions(self):
        data = np.genfromtxt(self.layout_file)
        N = int(sum(data[:-1]))
        xy_position = np.zeros((N, 2), dtype=np.float32)
        cr_position = np.zeros((N, 2), dtype=np.int32)
        ind_pos = 0
        individual = data[:-1]

        for ind in range(self.rows * self.cols):
            if individual[ind] == 1:
                r_i = np.floor(ind / self.cols)
                c_i = np.floor(ind - r_i * self.cols)
                cr_position[ind_pos, 0] = c_i
                cr_position[ind_pos, 1] = r_i
                xy_position[ind_pos, 0] = c_i * self.cell_width + 0.5 * self.cell_width
                xy_position[ind_pos, 1] = r_i * self.cell_height + 0.5 * self.cell_height
                ind_pos += 1

        return xy_position, cr_position

    def initialize_graph(self):
        N = self.xy_position.shape[0]

        for i in range(N):
            node_type = -1 if i in [7, 21, 36, 43] else 0
            self.G.add_node(i, tp=node_type)
            for j in range(N):
                if j != i:
                    self.G.add_edge(i, j)

        self.G = self.G.to_directed()
        self.pos = dict(enumerate(self.xy_position))

        for i in self.G.nodes():
            node_type = self.G.nodes[i]['tp']
            if node_type == 1:
                self.node_labels[i] = 'S' + str(i)
            elif node_type == 0:
                self.node_labels[i] = 'T' + str(i)
            elif node_type == -1:
                self.node_labels[i] = 'Sbs' + str(i)

    def create_optimization_model(self):
        self.model = gb.Model()
        self.flow_vars = self.model.addVars(self.G.edges(), lb=0.0, vtype=gb.GRB.CONTINUOUS, name='f')
        self.edge_vars = self.model.addVars(self.G.edges(), lb=0.0, ub=1.0, vtype=gb.GRB.BINARY, name='y')
        self.cable_vars = self.model.addVars(
            [(i, j, cable) for i, j in self.G.edges() for cable in range(len(self.cables))],
            lb=0.0, ub=1.0,
            obj=[cable['unit'] * self.distance(self.pos[i], self.pos[j]) for i, j in self.G.edges() for cable in
                 self.cables],
            vtype=gb.GRB.BINARY, name='x'
        )
        self.add_constraints()

    def add_constraints(self):
        self.model.update()

        for i, j in self.G.edges():
            self.model.addConstr(
                gb.quicksum(self.cable_vars[i, j, cable] for cable in range(len(self.cables))) == self.edge_vars[i, j])

        for h in self.G.nodes():
            node_type = self.G.nodes[h]['tp']
            if node_type == 0:
                self.model.addConstr(self.flow_vars.sum(h, '*') - self.flow_vars.sum('*', h) == 5)
            elif node_type == 1:
                self.model.addConstr(self.flow_vars.sum(h, '*') - self.flow_vars.sum('*', h) == 0)

        for i, j in self.G.edges():
            self.model.addConstr(
                gb.quicksum(self.cables[cable]["flow_max"] * self.cable_vars[i, j, cable] for cable in
                            range(len(self.cables))) >= self.flow_vars[i, j]
            )

        self.model.addConstrs((self.edge_vars.sum(h, '*') == 1 for h in self.G.nodes() if self.G.nodes[h]['tp'] == 0))
        self.model.addConstrs((self.edge_vars.sum(h, '*') == 0 for h in self.G.nodes() if self.G.nodes[h]['tp'] == -1))
        self.model.addConstrs((self.edge_vars.sum(h, '*') <= 1 for h in self.G.nodes() if self.G.nodes[h]['tp'] == 1))
        self.model.addConstrs((self.edge_vars.sum('*', h) <= 1 for h in self.G.nodes() if self.G.nodes[h]['tp'] == 1))
        self.model.addConstrs(
            (self.edge_vars.sum('*', h) <= self.max_cable_subst for h in self.G.nodes() if self.G.nodes[h]['tp'] == -1))

        collisions = self.find_collisions()
        for cross in collisions:
            self.model.addConstr(
                self.edge_vars[cross[0]] + self.edge_vars[cross[1]] + self.edge_vars[(cross[0][1], cross[0][0])] +
                self.edge_vars[(cross[1][1], cross[1][0])] <= 1
            )

    def find_collisions(self):
        collisions = []
        for i, j in self.G.edges():
            for k, h in self.G.edges():
                if k != i and h != j and k != j and h != i:
                    if self.intersect(self.pos[i], self.pos[j], self.pos[k], self.pos[h]):
                        collisions.append([(i, j), (k, h)])
        return collisions

    def distance(self, P1, P2):
        return math.sqrt((P1[0] - P2[0]) ** 2 + (P1[1] - P2[1]) ** 2)

    def ccw(self, A, B, C):
        return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0])

    def intersect(self, A, B, C, D):
        return self.ccw(A, C, D) != self.ccw(B, C, D) and self.ccw(A, B, C) != self.ccw(A, B, D)

    def optimize(self):
        self.model.update()
        self.model.optimize()
        return self.model
