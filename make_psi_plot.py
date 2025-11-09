import src.solver.define_equilibria as define_equilibria

import src.plotting.plot_sensors as plot_sensors


define_equilibria.define_multiple_equilibria(plot_quantities=True)

plot_sensors.plot_sensors_on_machine_geometry()