from AST2000SolarSystem import AST2000SolarSystem as ast
import matplotlib.pyp

user = 'markusbp'

seed = ast.get_seed(user)
solar_system = ast(seed)

x0 = solar_system.x0
y0 = solar_system.y0

plt.scatter(x0, y0, 'x', 'y', 'Initial Solar System Configuration')
