import streamlit as st
import sympy as sp
import numpy as np
import plotly.graph_objs as go

st.set_page_config(layout="wide")

# Define symbols
x = sp.symbols('x')
y = sp.Function('y')(x)
a0, a1 = sp.symbols('a0 a1')

# Function to compute coefficients
def compute_coefficients(N: int, P: sp.Expr, Q: sp.Expr, a0: float=a0, a1: float=a1) -> tuple[list, list]:
    """
    Computes power series coefficients for the differential equation y'' + P(x)y' + Q(x)y = 0.
    
    Args:
        N (int): Order of approximation.
        P (sp.Expr): Polynomial P(x).
        Q (sp.Expr): Polynomial Q(x).
        a0 (symbol): Initial coefficient.
        a1 (symbol): First-order coefficient.
    
    Returns:
        tuple: (list of computed coefficients, list of recurrence equations)
    """
    coefficients = [a0, a1]
    
    # Get the power series expansions for P(x) and Q(x)
    try:
        p_series = sp.series(P, x, 0, N).removeO().as_poly(x).all_coeffs()[::-1]
        q_series = sp.series(Q, x, 0, N).removeO().as_poly(x).all_coeffs()[::-1]
    except:
        raise ValueError("P(x) and Q(x) must be valid polynomials.")
    
    # Pad p_series and q_series to ensure they are at least length N
    p_series += [0] * (N - len(p_series))
    q_series += [0] * (N - len(q_series))
    
    # Symbols for new coefficients
    a = sp.symbols(f'a2:{N+2}')
    coefficients.extend(a)
    
    # Compute coefficients a_2 to a_N using the recurrence relation
    equations = []
    for n in range(N):
        lhs = (n + 2) * (n + 1) * coefficients[n + 2]
        rhs = -sum((k + 1) * p_series[n - k] * coefficients[k + 1] + q_series[n - k] * coefficients[k] for k in range(n + 1))
        equations.append(sp.Eq(lhs, rhs))
    
    # Solve the equations for a2, a3, ..., aN
    solutions = sp.solve(equations, a)
    
    # Update the coefficients with the solutions
    for sym, val in solutions.items():
        index = int(str(sym)[1:])  # get the numeric part after 'a'
        coefficients[index] = val
    
    return coefficients[:N + 2], equations

# Streamlit app
st.title("Power Series Method for Differential Equations")

st.latex(r"y'' + P(x) y' + Q(x) y = 0 \quad \text{(1)}")

# User input for polynomials P(x) and Q(x)
P_input = st.text_input("Enter polynomial P(x):", "x + 1")
Q_input = st.text_input("Enter polynomial Q(x):", "x")

# Convert input strings to SymPy expressions
try:
    P = sp.sympify(P_input)
except:
    P = None
try:
    Q = sp.sympify(Q_input)
except:
    Q = None 

# Check that P and Q are polynomials
if P is None or not P.is_polynomial(x):
    st.error("P(x) must be a polynomial. Please enter a valid polynomial.")
elif Q is None or not Q.is_polynomial(x):    
    st.error("Q(x) must be a polynomial. Please enter a valid polynomial.")
else:
    # Select N for the order of coefficients
    N = st.slider("Select the order N for coefficients calculation:", 2, 10, 4) - 1

    # Calculate coefficients and recurrence equations
    coeffs, equations = compute_coefficients(N, P, Q)

    # Display recurrence relations
    st.write("### Recursion Formula:")
    st.latex(r"(n+2)(n+1)a_{n+2} = - \sum_{k=0}^{n} \left[(k+1)p_{n-k} a_{k+1} + q_{n-k} a_{k}\right] \quad \text{(2)}")
    for eq in equations:
        st.latex(sp.latex(eq))

    # Display calculated coefficients
    st.write("### Coefficients:")
    for i, coeff in enumerate(coeffs):
        st.latex(f"a_{{{i}}} = {sp.latex(coeff)}")

    # Construct y as a power series
    y_series = sum(coeff * x**i for i, coeff in enumerate(coeffs))
    st.write("### Power Series Solution:")
    st.latex(f"y = {sp.latex(y_series)}")
   

    # Display general solution if possible
    st.write("### General Solution:")
    general_solution = sp.dsolve(sp.Derivative(y, x, x) + P * sp.Derivative(y, x) + Q * y, y)
    st.latex(sp.latex(general_solution))

    # Attempt to process general solution constants
    try:
        general_solution_rhs = general_solution.rhs.removeO()
        gen_constants = [symbol for symbol in general_solution_rhs.free_symbols if symbol != x]
        
        # Set up default values for constants (if available)
        gen_constant_values = {symbol: 1.0 for symbol in gen_constants}
        
        # Substitute constant values into the general solution
        general_solution_with_values = general_solution_rhs.subs(gen_constant_values)
        general_solution_func = sp.lambdify(x, general_solution_with_values, 'numpy')
        general_solution_available = True
    except:
        general_solution_available = False

    # Organize input fields into columns
    col1, col2 = st.columns(2)

    with col1:
        # Input fields for a0 and a1
        a0_val = st.number_input("Enter the value for a0:", value=1.0)
        a1_val = st.number_input("Enter the value for a1:", value=1.0)

    if general_solution_available:
        with col2:
            # Display input fields for each general constant (C1, C2, etc.)
            for const in gen_constants:
                gen_constant_values[const] = st.number_input(f"Enter the value for {const}:", value=1.0)

    # Generate x values and evaluate the general solution (if available)
    x_vals = np.linspace(-1, 1, 100)
    fig = go.Figure()

    # Plot partial sums of the series solution

    def get_series_approx(
        x_vals: np.ndarray, 
        coeffs: list[sp.Expr], 
        num_terms: int, 
        a0_val: float, 
        a1_val: float
    ) -> list[float]:
        """
        Computes partial sums of the power series solution.

        Args:
            x_vals (np.ndarray): Array of x-values where the series is evaluated.
            coeffs (List[Expr]): List of symbolic coefficients of the power series.
            num_terms (int): Number of terms to include in the series approximation.
            a0_val (float): Value of coefficient a0.
            a1_val (float): Value of coefficient a1.

        Returns:
            List[float]: List of evaluated y-values corresponding to x_vals.
        """
        y_vals = []
        for x_val in x_vals:
            y_val = sum(float(coeff.subs({a0: a0_val, a1: a1_val})) * (x_val ** i) for i, coeff in enumerate(coeffs[:num_terms]))
            y_vals.append(y_val)
        return y_vals

    for num_terms in range(1, len(coeffs) + 1):
        y_vals = get_series_approx(x_vals, coeffs, num_terms, a0_val, a1_val)
        fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode='lines', name=f'N={num_terms - 1}'))

    if general_solution_available:
        # Update general solution with user-defined constants and plot it
        try:
            updated_general_solution_with_values = general_solution_rhs.subs(gen_constant_values)
            general_solution_func = sp.lambdify(x, updated_general_solution_with_values, 'numpy')
            gen_solution_y_vals = general_solution_func(x_vals)

            if np.all(gen_solution_y_vals == 0):
                gen_solution_y_vals = np.zeros_like(x_vals)  
            
            fig.add_trace(go.Scatter(x=x_vals, y=gen_solution_y_vals, mode='lines', name="General Solution",
                                    line=dict(color='black', dash='dash')))
        except Exception as e:
            st.error(f"ERROR: Invalid expression for the general solution {general_solution_with_values}: {e}")
    
    # Plot settings
    fig.update_layout(
        height=700,
        xaxis_title="x",
        yaxis_title="y(x)",
        legend=dict(font=dict(size=18)),
        annotations=[
            dict(
                xref="paper", yref="paper",
                x=0.5, y=1.15,
                text=f"a0 = {a0_val}, a1 = {a1_val}",
                showarrow=False,
                font=dict(size=18)
            ),
            dict(
                xref="paper", yref="paper",
                x=0.5, y=1.10,
                text=", ".join([f"{const} = {gen_constant_values[const]}" for const in gen_constants]),
                showarrow=False,
                font=dict(size=18)
            ) if general_solution_available else {}
        ]
    )

    fig.update_xaxes(zeroline=True, zerolinewidth=2, showgrid=True)
    fig.update_yaxes(zeroline=True, zerolinewidth=2, showgrid=True)

    # Show the plot    
    st.plotly_chart(fig, use_container_width=True)
