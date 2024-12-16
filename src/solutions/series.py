import sympy as sp

def get_derivatives(order, n):

    x = sp.symbols('x')
    
    a = sp.symbols(f'a0:{n}')  # Coefficients a0, a1, ..., an

    # Taylor series: y(x) = a0 + a1*x + a2*x**2 + ... + an*x**n
    series = sum(a[n] * (x ** n) for n in range(n))
    print("series:", series) 

    # Compute derivatives of the Taylor series up to the required order
    derivatives = [series]
    for n in range(1, order + 1):
        derivatives.append(sp.diff(derivatives[-1], x))

    return derivatives
