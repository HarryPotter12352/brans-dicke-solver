import sympy as sp

# Define some useful/commonly used symbols. Note: Do NOT enter a symbol as a coordinate.
a, k, rho, p = sp.symbols("a k rho p")

# Collecting all inputs
def input_coordinates():
    coords = input("Enter the coordinates (comma-separated): ").split(",")
    return sp.symbols(','.join(coords))

def input_tensor(type: str, size: int):
    print(f"Enter the elements of the {size}x{size} {type} tensor")
    tensor_elements = []
    for i in range(size):
        row = input(f"Row {i+1} (comma-separated): ").split(",")
        tensor_elements.append([sp. sympify(element.replace("^", "**")) for element in row])
    return sp.Matrix(tensor_elements)

def input_scalar_field():
    return sp.sympify(input("Enter the scalar field: ").replace("^", "**"))

def input_potential():
    return sp.sympify(input("Enter the potential: ").replace("^", "**"))


coordinates = input_coordinates()
size = len(coordinates)

omega = int(input("Enter the coupling factor (integer): "))


metric = input_tensor("metric", size)
energy_momentum_tensor = input_tensor("energy-momentum", size)
potential = input_potential()
derivative_of_potential = potential.diff("phi")

phi = input_scalar_field()
derivative_of_potential = derivative_of_potential.subs(sp.Symbol("phi"), phi)
derivatives_of_scalar_field = [phi.diff(coord) for coord in coordinates]

inverse = metric.inv()
determinant = metric.det()
# Compute the kinetic term g^{\alpha\beta}\partial_{\alpha}\phi\partial_{\beta}\phi
kinetic_term = 0
for mu in range(size):
    for nu in range(size):
        kinetic_term += inverse[mu, nu] * derivatives_of_scalar_field[mu] * derivatives_of_scalar_field[nu]

# Compute the trace of the energy momentum tensor, given by T^{\alpha}_{\alpha} = g^{\alpha\beta}T_{\alpha\beta}
trace_of_energy_momentum_tensor = 0
for mu in range(size):
    for nu in range(size):
        trace_of_energy_momentum_tensor += inverse[mu, nu] * energy_momentum_tensor[mu, nu]

# Compute the Laplace-Beltrami operator, given by \Box\phi = \frac{8\piT + 2V - \phi V'(\phi)}{3 + 2\omega}
laplace_beltrami = ((8 * sp.pi * trace_of_energy_momentum_tensor) + 
                    (2 * potential - phi * derivative_of_potential)) / (3 + 2 * omega)

# Finally, put everything together to compute the Einstein tensor
einstein_tensor = sp.Matrix.zeros(size, size)
for mu in range(size):
    for nu in range(size):
        einstein_tensor[mu, nu] = ((8 * sp.pi / phi) * trace_of_energy_momentum_tensor + 
                                   (omega / phi**2) * (derivatives_of_scalar_field[mu] * derivatives_of_scalar_field[nu] - 
                                   0.5 * metric[mu, nu] * kinetic_term) + 
                                   (1 / phi) * (phi.diff(coordinates[mu], coordinates[nu]) - 
                                   metric[mu, nu] * laplace_beltrami) - 
                                   0.5 * metric[mu, nu] * potential / (phi))

print(f"Einstein Tensor:\n{einstein_tensor}")