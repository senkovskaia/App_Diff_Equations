import streamlit as st
import sympy as sp
from sympy import symbols, Eq, solve, sympify
import plotly.graph_objects as go
import numpy as np

st.set_page_config(layout="wide")

# Define x symbol globally
x = symbols('x')

def compute_a_coefficients(p_dict: dict[int, sp.Expr], 
                           q_dict: dict[int, sp.Expr], 
                           r: sp.Symbol, 
                           n_max: int, 
                           a0_val: sp.Expr
                           ) -> dict[int, sp.Expr]:
    """
    Computes the coefficients a_n recursively up to n_max.
    Args:
        p_dict (dict): Dictionary of coefficients for P(x).
        q_dict (dict): Dictionary of coefficients for Q(x).
        r (sympy.Expr): Indicial root.
        n_max (int): Maximum order for computation.
        a0_val (float): Initial coefficient a0.

    Returns:
        dict: Computed coefficients a_n.
    """
    a = symbols('a_0:%d' % (n_max + 1))
    
    # Base case: a_0
    a_values = {0: a0_val}
    
    # Compute a_n recursively
    for n in range(1, n_max + 1):
        
        first_tirm = a[n] * ((r + n) * (r + n - 1) + (r + n) * p_dict[0] + q_dict[0])

        summation_term = sum(a_values[k] * ((r + k) * p_dict[n - k] + q_dict[n - k]) for k in range(n))

        recurrence_eq = Eq(first_tirm + summation_term, 0)
        
        solution = solve(recurrence_eq, a[n])
        
        if not solution:
            a_n_value = 0
        else:
            a_n_value = solution[0]
        
        a_values[n] = a_n_value
    
    return a_values
    

def convert_to_integer_keys(input_dict: dict[sp.Symbol, sp.Expr], n: int) -> dict[int, sp.Expr]:
    """
    Converts a dictionary with keys as powers of x (e.g., x, x^2) to one with integer keys.

    Args:
        input_dict (dict): SymPy dictionary of coefficients.
        n (int): Maximum power to consider.

    Returns:
        dict: Dictionary with integer keys representing polynomial degrees.
    """
    x = symbols('x')
    output_dict = {}
    
    # Fill the dictionary for all powers up to max_power
    for power in range(n+1): 
        output_dict[power] = 0  # Default value for missing powers

    for key, value in input_dict.items():
        if key == 1:  # Constant term (x^0)
            output_dict[0] = value
        else:
            # Extract the power of x from the symbolic key
            power = key.as_poly(x).degree() if key.has(x) else 0
            output_dict[power] = value
    
    return output_dict


# Title and Explanation
st.title("Frobenius Method for Differential Equations")

# Display Formulas
st.latex(r"y'' + P(x)y' + Q(x)y = 0 \quad \text{(1)}")
st.latex(r"y = x^r(a_0 + a_1x + a_2x^2 + \dots) \quad \text{(2)}")
st.latex(r"xP(x) = \sum_{n=0}^\infty p_n x^n, \quad x^2Q(x) = \sum_{n=0}^\infty q_n x^n \quad \text{(3)}")

# Input Fields for P(x) and Q(x)
col1, col2 = st.columns([0.25, 2])  # First column narrower, second column wider

# First row for xP(x)
with col1:
    st.markdown("", unsafe_allow_html=True)  # Adds a small vertical space
    st.latex(r"xP(x) =")
with col2:
    xP_input = st.text_input("Enter polynomial xP(x):", "x + 2")

# Second row for x^2Q(x)
with col1:
    st.markdown("######", unsafe_allow_html=True)  # Adds consistent spacing
    st.latex(r"x^2Q(x) =")
with col2:
    x2Q_input = st.text_input("Enter polynomial x^2Q(x):", "x")

# Convert input to SymPy expressions
try:
    P = sympify(xP_input)
    Q = sp.sympify(x2Q_input)
except:
    P = None
    Q = None

# Validate polynomials
if P is None or not P.is_polynomial(x):
    st.error("xP(x) must be a polynomial. Please enter a valid polynomial.")
elif Q is None or not Q.is_polynomial(x):
    st.error("x^2Q(x) must be a polynomial. Please enter a valid polynomial.")
else:
    # Select N for the order of coefficients
    N = st.slider("Select the maximum order N for coefficients calculation:", 2, 10, 4) - 1
    # Display the entered equation
    st.write("### Differential equation:")
    st.latex(r"y'' + \frac{%s}{x}y' + \frac{%s}{x^2}y = 0" % (sp.latex(P), sp.latex(Q)))
    # Extract p0 and q0
    P_coeff_dict = P.as_coefficients_dict()
    Q_coeff_dict = Q.as_coefficients_dict()

    p0 = P_coeff_dict.get(1, 0) 
    q0 = Q_coeff_dict.get(1, 0)

    p_dict = convert_to_integer_keys(P_coeff_dict, N + 1)
    q_dict = convert_to_integer_keys(Q_coeff_dict, N + 1)

    # Solve the indicial equation
    r = sp.symbols('r')
    indicial_eq = r * (r - 1) + p0 * r + q0
    roots = sp.solve(indicial_eq, r)

    if not all(root.is_real for root in roots):
        #st.error(f"Warning: complex numbers (!)")
        pass
    else:
        roots = sorted(roots, reverse=True) #  [0, -1]

    st.write("### Recursion formula:")
    st.latex(r"""
            a_n \left[(r+n)(r+n-1) + (r+n)p_0 + q_0 \right] + 
            \sum_{k=0}^{n-1} a_k \left[(r+k)p_{n-k} + q_{n-k} \right] = 0 
            \quad \text{(4)}
            """)

    #st.write("The first equation (a0 ≠ 0):")
    st.latex(r"r(r-1) + rp_0 + q_0 = 0 \quad \text{(5)}")
    st.latex(r"p_0 = %s, \quad q_0 = %s" % (sp.latex(p0), sp.latex(q0)))
    st.latex(sp.latex(indicial_eq) + " = 0")
    # Display the roots
    r1 = roots[0]
    if len(roots) == 1:
        r2 = r1
    elif len(roots) == 2:
        r2 = roots[1]
    else:
        st.error("There should be only 1 or 2 roots")

    r_diff = r1 - r2

    st.latex(r"r_1 = %s, \quad r_2 = %s" % (sp.latex(r1), sp.latex(r2)))
    st.latex(r"r_1 - r_2 = %s" % (sp.latex(r_diff)))

    case_num = 0
    # value is a whole number
    if r_diff == 0:
        case_num = 3  
        st.write("### Case 3:")
        st.latex(r"\text{Since } r1 = r2 \text{, there is only one Frobenius solution:}")
        st.latex(r"y_1(x) = x^{r_1}\sum_{n=0}^\infty a^{(1)}_n x^n")
        st.latex(r"\text{The second independent solution has the form:}")
        st.latex(r"y_2(x) = y_1(x)log(x) + x^{r_2}\sum_{n=0}^\infty a^{(2)}_n x^n")

    elif r_diff % 1 == 0:
        case_num = 2
        st.write("### Case 2:")
        st.latex(r"\text{Since } r1 - r2 \in \mathbb{N} \text{, there is only one Frobenius solution:}")
        st.latex(r"y_1(x) = x^{r_1}\sum_{n=0}^\infty a^{(1)}_n x^n")
        st.latex(r"\text{The second independent solution has the form:}")
        st.latex(r"y_2(x) = C \cdot y_1(x)log(x) + x^{r_2}\sum_{n=0}^\infty a^{(2)}_n x^n")

    elif r_diff % 1 != 0:
        case_num = 1
        st.write("### Case 1:")
        st.latex(r"\text{Since } r1 - r2 \not\in \mathbb{Z} \text{, there are two independent solutions in Frobenius form:}")
        st.latex(r"y_1(x) = x^{r_1}\sum_{n=0}^\infty a^{(1)}_n x^n")
        st.latex(r"y_2(x) = x^{r_2}\sum_{n=0}^\infty a^{(2)}_n x^n")
    
    else:
        st.error("Something wrong!")

    st.latex(r'\text{General solution:}')
    st.latex('y(x) = C1*y_1(x) + C2*y_2(x)')
    
    a0_val = symbols('a_0') 

    a_values_1 = compute_a_coefficients(p_dict, q_dict, r1, N + 1, a0_val)
    a_values_2 = compute_a_coefficients(p_dict, q_dict, r2, N + 1, a0_val)


    st.write("### Coefficients:")
    
    col3, col4 = st.columns(2)

    with col3:
        for i, coeff in a_values_1.items():
            st.latex(r"a^{(1)}" +f"_{{{i}}} = {sp.latex(coeff)}")

    with col4:
        for i, coeff in a_values_2.items():
            st.latex(r"a^{(2)}" +f"_{{{i}}} = {sp.latex(coeff)}")


    if not all(root.is_real for root in roots):
        st.error("The plot cannot be generated for complex roots")
    else:

        # Add a reset button
        if "a0" not in st.session_state:
            st.session_state.a0 = 1.0
        if "C1" not in st.session_state:
            st.session_state.C1 = 1.0
        if "C2" not in st.session_state:
            st.session_state.C2 = 1.0
        if "C" not in st.session_state and case_num == 2:
            st.session_state.C = 1.0

        
        if st.button("Reset to Default"):
            st.session_state.a0 = 1.0
            st.session_state.C1 = 1.0
            st.session_state.C2 = 1.0
            if case_num == 2:
                st.session_state.C = 1.0

        # Sliders linked to session state
        a0 = st.slider("Select the value of a_0:", 1.0, 10.0, key="a0")
        C1_val = st.slider("Select the value of C1:", -10.0, 10.0, key="C1")
        C2_val = st.slider("Select the value of C2:", -10.0, 10.0, key="C2")

        if case_num == 2:
            C = st.slider("Select the value of C:", -10.0, 10.0, key="C")
        else:
            C = None 

        # Generate x values
        x_vals = np.linspace(0.1, 5, 500)
        x_vals = x_vals[x_vals > 0]  # Avoid x = 0 due to log(x)

        def compute_y_symbolic(a_values: dict[int, sp.Expr], r: sp.Symbol) -> sp.Expr:
            """
            Compute the Frobenius series solution y(x) for given coefficients.
            """
            x = symbols('x')
            y = 0  # Initialize y as a symbolic expression

            for n, coeff in a_values.items():
                # Build the Frobenius series term by term
                y += coeff * x**(n + r)
            
            return y

        y1_symbolic = compute_y_symbolic(a_values_1, r1)

        if case_num == 1:
            # Case 1: Simple Frobenius form
            y2_symbolic = compute_y_symbolic(a_values_2, r2)
        elif case_num == 2:
            # Case 3: Includes logarithmic term with constant C
            x = symbols('x')
            log_term = C * y1_symbolic * sp.log(x)
            frobenius_term = compute_y_symbolic(a_values_2, r2)
            y2_symbolic = log_term + frobenius_term
        elif case_num == 3:
            # Case 2: Includes logarithmic term with constant C1
            x = symbols('x')
            log_term = C1_val * y1_symbolic * sp.log(x)
            frobenius_term = compute_y_symbolic(a_values_2, r2)
            y2_symbolic = log_term + frobenius_term
        else:
            st.error("Unexpected case number!")

        # Convert symbolic expressions to numerical functions
        y1_func = sp.lambdify(symbols('x'), y1_symbolic.subs({symbols('a_0'): a0}), modules="numpy")
        y2_func = sp.lambdify(symbols('x'), y2_symbolic.subs({symbols('a_0'): a0}), modules="numpy")

        # Evaluate the numerical functions
        y1_vals = y1_func(x_vals)
        y2_vals = y2_func(x_vals)

        # Combine y1 and y2 with C1 and C2
        y_combined_vals = C1_val * y1_vals + C2_val * y2_vals

        # Plotting using Plotly
        fig = go.Figure()

        # Add traces for y1, y2, and y_combined
        fig.add_trace(go.Scatter(x=x_vals, y=y1_vals, mode='lines', name='y1'))
        fig.add_trace(go.Scatter(x=x_vals, y=y2_vals, mode='lines', name='y2'))
        fig.add_trace(go.Scatter(x=x_vals, y=y_combined_vals, mode='lines', name='y = C1*y1 + C2*y2',
                                 line=dict(color='black', dash='dash')))

        # Layout configuration
        fig.update_layout(
            xaxis_title="x",
            yaxis_title="y(x)",
            template="plotly_white",
            height=700,
            legend=dict(font=dict(size=18)),
            annotations=[
            dict(
                xref="paper", yref="paper",
                x=0.5, y=1.15,
                text=f"a0 = {a0}, C1 = {C1_val}, C2 = {C2_val}" + (f", C = {C}" if C else ""),
                showarrow=False,
                font=dict(size=18)
            )
        ]
        )
        fig.update_xaxes(zeroline=True, zerolinewidth=2, showgrid=True)
        fig.update_yaxes(zeroline=True, zerolinewidth=2, showgrid=True)

        st.plotly_chart(fig)