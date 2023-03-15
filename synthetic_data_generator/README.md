## MiADE Synthetic Data Generator and Trainer

A helper streamlit app for MiADE synthetic data generation

To run locally:
```bash
pip install -r requirements.txt
streamlit run app.py
```
In Docker:
```bash
 docker build -t streamlit .
 docker run --name miade-data-app -p 8501:8501 streamlit
```