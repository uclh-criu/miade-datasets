import pandas as pd
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder, ColumnsAutoSizeMode
from st_aggrid.shared import GridUpdateMode

data = pd.read_csv(
    "./samples/problems_synthetic_train_example.csv"
)

st.set_page_config(
    layout="wide", page_icon="🖱️", page_title="Interactive table app"
)
st.title("🖱️ MiADE MetaCAT training data")
st.write(
    """Interact with MiADE's train dataset"""
)

MIN_HEIGHT = 27
MAX_HEIGHT = 800
ROW_HEIGHT = 35


def aggrid_interactive_table(df: pd.DataFrame):
    """Creates an st-aggrid interactive table based on a dataframe.

    Args:
        df (pd.DataFrame]): Source dataframe

    Returns:
        dict: The selected row
    """
    options = GridOptionsBuilder.from_dataframe(
        df, enableRowGroup=True, enableValue=True, enablePivot=True, min_column_width=100
    )
    options.configure_selection(selection_mode="multiple", use_checkbox=True)
    options.configure_side_bar()

    options.configure_selection("single")
    selection = AgGrid(
        df,
        height=min(MIN_HEIGHT + len(df) * ROW_HEIGHT, MAX_HEIGHT),
        columns_auto_size_mode=ColumnsAutoSizeMode.FIT_ALL_COLUMNS_TO_VIEW,
        fit_columns_on_grid_load=True,
        gridOptions=options.build(),
        theme="streamlit",
        update_mode=GridUpdateMode.MODEL_CHANGED,
        allow_unsafe_jscode=True,
    )

    return selection


selection = aggrid_interactive_table(df=data)

if selection:
    st.write("You selected:")
    st.json(selection["selected_rows"])
