import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

def preprocess_input(equation_str: str) -> str:
    """
    Preprocess input string by replacing 1st and 2nd derivatives (y' and y'') with SymPy diff expressions.
    Handles both first and second derivatives for symbolic computation.
    """
    equation_str = equation_str.replace("y''", "Derivative(y, x, 2)") 
    equation_str = equation_str.replace("y'", "Derivative(y, x, 1)")
    return equation_str


def string_to_sympy_eq(equation_str: str) -> sp.Expr:
    """
    Convert a user-inputted equation string into a SymPy equation.
    """
    x = sp.symbols('x')
    y = sp.Function('y')(x)
    equation_str_proc = preprocess_input(equation_str)
    equation = parse_expr(equation_str_proc, transformations="all", local_dict={'x': x, 'y': y})
    return equation