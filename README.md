# 🧬 SeqAlign Pro

**SeqAlign Pro** is a bioinformatics web application that visualizes the **Needleman-Wunsch algorithm** for global sequence alignment. It provides interactive dot plots, dynamic scoring matrices, and optimal path traceback, making it an excellent tool for educational and research purposes.

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://seq-align-pro.streamlit.app/)

## 📌 Overview
**SeqAlign Pro** is designed to simplify pairwise sequence alignment. Unlike standard command-line tools, it offers a visual dashboard to understand *how* sequences align.

It handles **FASTA file uploads** and performs **Global Alignment (Needleman-Wunsch)**, allowing users to identify mutations, gaps, and conserved regions between two DNA sequences efficiently.

## 🚀 Key Features
* **Global Alignment:** Implements the classic Needleman-Wunsch algorithm with customizable scoring (Match +1, Mismatch -1, Gap -2).
* **Dot Plots:** Visualizes sequence similarity using **Altair** heatmaps within scrollable containers.
* **Traceback Reconstruction:** Backtracks through the dynamic programming matrix to find the optimal alignment path.
* **Bulk Data Loading:** Supports FASTA file uploads to load sequences instantly.
* **Complexity Management:** Includes degradation logic to safely handle larger datasets without crashing the browser.

## 🛠️ Tech Stack
* **Language:** Python 3.9+
* **Framework:** Streamlit (Web UI)
* **Computation:** NumPy (Vectorized operations)
* **Data Manipulation:** Pandas
* **Visualization:** Altair

## 💻 How to Run Locally
If you want to run this tool on your own machine:

1. **Clone the repository**
    ```bash
    git clone [https://github.com/rao-monika-yadav/sequence-alignment.git](https://github.com/rao-monika-yadav/sequence-alignment.git)
    ```
2. **Install dependencies**
    ```bash
    pip install streamlit pandas numpy altair
    ```
3. **Run the app**
    ```bash
    streamlit run app.py
    ```

## 👨‍🔬 Author
**Monika Yadav**
* MSc Bioinformatics Candidate
* BSc (Hons.)/MSc Biotechnology

[LinkedIn](will add soon) | [GitHub](https://github.com/rao-monika-yadav)

---

*Disclaimer: This tool is intended for research and educational purposes only.*
