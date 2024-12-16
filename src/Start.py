import streamlit as st

st.set_page_config(
    page_title="calculator",
    page_icon="👋",
)

st.write("## Welcome to the Differential Equations App! 👋")

st.sidebar.success("Select a page above.")

st.markdown(
    """
    This application is designed to demonstrate solutions to differential equations using power series.
    
    **👈 Select a page from the sidebar** 
"""
)