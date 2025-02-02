import sympy as sp

def get_derivatives(order: int, n: int) -> list[sp.Expr]:
    """
    Computes derivatives of a Taylor series expansion up to the given order.
    
    Args:
        order (int): Number of derivatives required.
        n (int): Number of terms in the Taylor series.
    
    Returns:
        list: A list of SymPy expressions representing the derivatives.
    """
    x = sp.symbols('x')
    a = sp.symbols(f'a0:{n}') # Coefficients a0, a1, ..., an

    # Taylor series: y(x) = a0 + a1*x + a2*x**2 + ... + an*x**n
    series = sum(a[n] * (x ** n) for n in range(n))

    # Compute derivatives of the Taylor series up to the required order
    derivatives = [series]
    for n in range(1, order + 1):
        derivatives.append(sp.diff(derivatives[-1], x))

    return derivatives
