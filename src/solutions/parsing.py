import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

def preprocess_input(equation_str):
    # Preprocessing input to replace y' and y'' with diff expressions
    equation_str = equation_str.replace("y''", "Derivative(y, x, 2)")      # Replace y' with first derivative
    equation_str = equation_str.replace("y'", "Derivative(y, x, 1)")  # Replace y'' with second derivative
    return equation_str

# Function to convert user input string to SymPy equation
def string_to_sympy_eq(equation_str):
    x = sp.symbols('x')
    y = sp.Function('y')(x)
    equation_str_proc = preprocess_input(equation_str)
    equation = parse_expr(equation_str_proc, transformations="all", local_dict={'x': x, 'y': y})
    return equation