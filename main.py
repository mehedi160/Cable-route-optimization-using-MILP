from cable_routing_MILP import CableRoutingMILP
from plot_route import plot_solution

# Initialize the optimizer
layout_file = 'best_layouts.dat'
wind_farm_optimizer = WindFarmOptimizer(layout_file)

# Run the optimization
model = wind_farm_optimizer.optimize()

# Plot the solution
plot_solution(
    wind_farm_optimizer.G,
    wind_farm_optimizer.pos,
    wind_farm_optimizer.cable_vars,
    wind_farm_optimizer.cables,
    wind_farm_optimizer.xy_position,
    wind_farm_optimizer.node_labels
)
