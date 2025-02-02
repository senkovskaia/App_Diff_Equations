## Differential Equations App

This a **Streamlit** application that is part of the **bachelor's thesis of Klavdiia Senkovskaia**, designed to interactively demonstrate solutions to differential equations using power series. It focuses on providing an interactive approach to solving differential equations using:

1. **Power Series Method**

2. **Frobenius Method**

It provides interactive inputs, symbolic computations using **SymPy**, and visual plots using **Plotly**.

## Features
- Compute coefficients for power series solutions using symbolic differentiation.
- Solve linear differential equations with variable coefficients.
- Recurrence relations for finding coefficients of power series expansions.
- Interactive visualization of solutions using numerical approximation.
- Handles different cases of the Frobenius method, including singular points.

## Hosted Application

This application is hosted on **Streamlit Community Cloud** and can be accessed here:

🔗 [https://app-diff-equations.streamlit.app](https://app-diff-equations.streamlit.app)

## Installation

### Requirements
This project requires **Python 3.12+** and uses **Poetry** for dependency management.
### 1. Clone the Repository
```sh
git clone https://github.com/senkovskaia/App_Diff_Equations
cd App_Diff_Equations

```

### 2. Install Dependencies
Make sure you have **Python 3.12+** installed. Then run:
```sh
poetry install
poetry shell
```

### 3. Run the Application
```sh
streamlit run src/Start.py
```

## File Structure
```
App_Diff_Equations/
├── pyproject.toml       # Python package configuration
├── README.md      
├── src/
│   ├── Start.py         # Main Streamlit app
│   ├── solutions/
│   │   ├── parsing.py   # Expression parsing utilities
│   │   ├── series.py    # Power series computations
│   ├── pages/
│   │   ├── 1_PowerSeries.py   # Power series method implementation
│   │   ├── 2_Frobenius.py     # Frobenius method implementation
```

---
### Contributors
Developed by **Klavdiia Senkovskaia**. Feel free to contribute!

For issues & feedback, open an **issue** in the repository.
