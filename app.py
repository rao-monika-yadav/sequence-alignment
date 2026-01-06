import streamlit as st
import numpy as np
import pandas as pd
import altair as alt
from io import StringIO

# --- PART 1: BACKEND LOGIC ---

def create_dotplot(seq1, seq2):
    s1 = np.array(list(seq1.upper()))
    s2 = np.array(list(seq2.upper()))
    return (s1 == s2[:, None]).astype(int)

def needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-2):
    seq1, seq2 = seq1.upper(), seq2.upper()
    n, m = len(seq1), len(seq2)
    score_matrix = np.zeros((m + 1, n + 1))
    
    for i in range(m + 1): score_matrix[i][0] = i * gap_penalty
    for j in range(n + 1): score_matrix[0][j] = j * gap_penalty
        
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score = match_score if seq1[j-1] == seq2[i-1] else mismatch_penalty
            score_matrix[i][j] = max(
                score_matrix[i-1][j-1] + score,
                score_matrix[i-1][j] + gap_penalty,
                score_matrix[i][j-1] + gap_penalty
            )
    return score_matrix

def run_traceback(seq1, seq2, score_matrix, match_score=1, mismatch_penalty=-1, gap_penalty=-2):
    seq1, seq2 = seq1.upper(), seq2.upper()
    i, j = len(seq2), len(seq1)
    align1, align2 = "", ""
    while i > 0 or j > 0:
        curr = score_matrix[i][j]
        if i > 0 and j > 0:
            penalty = match_score if (seq1[j-1] == seq2[i-1]) else mismatch_penalty
            if curr == score_matrix[i-1][j-1] + penalty:
                align1 += seq1[j-1]; align2 += seq2[i-1]; i -= 1; j -= 1; continue
        if i > 0 and curr == score_matrix[i-1][j] + gap_penalty:
            align1 += "-"; align2 += seq2[i-1]; i -= 1; continue
        if j > 0 and curr == score_matrix[i][j-1] + gap_penalty:
            align1 += seq1[j-1]; align2 += "-"; j -= 1; continue
    return align1[::-1], align2[::-1]

def plot_dotplot(matrix, seq1, seq2):
    """
    Visualizes the Dot Plot (Static View).
    """
    df = pd.DataFrame(matrix).stack().reset_index()
    df.columns = ['y', 'x', 'match']
    matches = df[df['match'] == 1].copy()
    matches['base'] = matches['x'].apply(lambda x: list(seq1.upper())[x])
    
    # Color Scale
    color_scale = alt.Color(
        'base:N', 
        scale=alt.Scale(
            domain=['A', 'T', 'G', 'C', 'N'], 
            range=['#50C878', '#E74C3C', '#F39C12', '#3498DB', '#95A5A6']
        ),
        legend=None 
    )
    
    chart = alt.Chart(matches).mark_rect().encode(
        x=alt.X('x:O', axis=alt.Axis(labelAngle=0, orient='top'), title='Sequence 1 (Horizontal)'),
        y=alt.Y('y:O', title='Sequence 2 (Vertical)'),
        color=color_scale,
        tooltip=['x', 'y', 'base']
    ).properties(width=len(seq1)*40, height=len(seq2)*40)
    
    return chart

# --- PART 2: UI CONFIGURATION ---

st.set_page_config(page_title="SeqAlign Pro", page_icon="🧬", layout="wide", initial_sidebar_state="expanded")

# Initialize Session States
if 'seq1' not in st.session_state: st.session_state['seq1'] = "GATTACA"
if 'seq2' not in st.session_state: st.session_state['seq2'] = "GCATGC"

# --- SIDEBAR ---
with st.sidebar:
    st.markdown("""<div style="text-align: center;"><div style="display: inline-flex; justify-content: center; align-items: center; width: 60px; height: 60px; background-color: #0068c9; border-radius: 50%; box-shadow: 2px 2px 5px rgba(0,0,0,0.2); margin-bottom: 15px;"><span style="color: white; font-size: 24px; font-weight: bold; font-family: 'Lucida Handwriting', cursive;">MY</span></div></div>""", unsafe_allow_html=True)
    st.title("SeqAlign Pro")
    st.caption("v1.0.0 | Bio-Tools Suite")
    st.markdown("---")
    
    st.write("### 🛠️ Key Features")
    st.markdown("""
    - **Global Alignment:** Uses the Needleman-Wunsch algorithm for optimal global matching.
    - **Dot Plots:** Interactive heatmaps to visualize sequence similarity and repeats.
    - **Traceback:** Reconstructs the exact alignment path (including gaps/mismatches).
    """)
    
    st.markdown("---")
    st.write("### 👨‍🔬 Author"); st.markdown("**Monika Yadav**"); st.caption("MSc Bioinformatics Candidate"); st.caption("BSc (Hons.)/MSc Biotechnology")
    st.markdown("""<div style="display: flex; gap: 10px;"><a href="#will add soon" target="_blank">LinkedIn</a> • <a href="https://github.com/rao-monika-yadav" target="_blank">GitHub</a> • <a href="# will add soon">Email</a></div>""", unsafe_allow_html=True)
    st.markdown("---")
    st.warning("For Educational and Research Use Only")

# --- MAIN WORKSPACE ---
st.title("🧬 SeqAlign Pro")

with st.expander("📖 How to use this tool"):
    st.markdown("1. **Load Data:** Upload FASTA or type manually.\n2. **Analyze:** Click the button.\n3. **Explore:** Use the Dot Plot (Tab 1) or Alignment (Tab 2).")

# --- A. DATA LOADING ---
with st.expander("📂 Load Data from FASTA File"):
    uploaded_file = st.file_uploader("Upload .fasta or .txt file", type=["fasta", "txt"])
    if uploaded_file:
        stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
        seqs = []
        curr = ""
        for line in stringio:
            line = line.strip()
            if line.startswith(">"):
                if curr: seqs.append(curr)
                curr = ""
            else: curr += line
        if curr: seqs.append(curr)
        
        if len(seqs) >= 2:
            st.session_state['seq1'] = seqs[0]
            st.session_state['seq2'] = seqs[1]
            if len(seqs) > 2: st.warning(f"⚠️ File contained {len(seqs)} sequences. Loaded first 2.")
            else: st.success("✅ Loaded 2 sequences successfully!")
        elif len(seqs) == 1:
            st.session_state['seq1'] = seqs[0]
            st.info("ℹ️ Loaded 1 sequence. Enter second manually.")

# --- B. INPUTS ---
st.markdown("### Input Sequences")
c1, c2 = st.columns(2)
with c1: st.text_area("Sequence 1 (Horizontal)", height=100, key="seq1")
with c2: st.text_area("Sequence 2 (Vertical)", height=100, key="seq2")

run_btn = st.button("Analyze Alignment", type="primary", use_container_width=True)

# --- C. RESULTS ---
if run_btn:
    s1 = st.session_state['seq1'].strip().upper()
    s2 = st.session_state['seq2'].strip().upper()
    
    MAX_CALC_LENGTH = 1500    # to prevent crashes
    MAX_VIZ_LENGTH = 300
    
    if not s1 or not s2:
        st.error("Please provide both sequences.")
    elif len(s1) > MAX_CALC_LENGTH or len(s2) > MAX_CALC_LENGTH:
        st.error(f"⚠️ Input too large! Max limit is {MAX_CALC_LENGTH} bases.")
    else:
        st.divider()
        is_large_data = len(s1) > MAX_VIZ_LENGTH or len(s2) > MAX_VIZ_LENGTH
        
        if is_large_data:
            st.warning(f"⚠️ Large sequences detected (> {MAX_VIZ_LENGTH} bp). Visualizations hidden.")
            t2, = st.tabs(["📊 Alignment Result"])
            t1 = None
        else:
            t1, t2 = st.tabs(["🔍 Dot Plot", "📊 Alignment Result"])
        
        # --- TAB 1: DOT PLOT ---
        if t1 and not is_large_data:
            with t1:
                st.subheader("Sequence Dot Plot")
                
                # EXTERNAL LEGEND (Always visible)
                st.markdown("""
                <div style="display: flex; align-items: center; gap: 15px; margin-bottom: 10px;">
                    <span style='font-weight:bold; margin-right:10px;'>Key:</span>
                    <div><span style='color:#50C878; font-size:1.2em;'>●</span> A</div>
                    <div><span style='color:#E74C3C; font-size:1.2em;'>●</span> T</div>
                    <div><span style='color:#F39C12; font-size:1.2em;'>●</span> G</div>
                    <div><span style='color:#3498DB; font-size:1.2em;'>●</span> C</div>
                    <div><span style='color:#95A5A6; font-size:1.2em;'>●</span> N</div>
                </div>
                """, unsafe_allow_html=True)

              #  st.caption("📜 Use the container scrollbars to view large sequences.")
                with st.container(height=500, border=True):
                    dotplot_chart = plot_dotplot(create_dotplot(s1, s2), s1, s2)
                    st.altair_chart(dotplot_chart, use_container_width=False)
            
        # --- TAB 2: ALIGNMENT ---
        with t2:
            st.subheader("Global Alignment Analysis")
            with st.spinner("Aligning sequences..."):
                mat = needleman_wunsch(s1, s2)
                res1, res2 = run_traceback(s1, s2, mat)
            
            if not is_large_data:
                with st.expander("View Scoring Matrix", expanded=True):
                    cols = [f"{c} ({i})" for i, c in enumerate(["-"]+list(s1))]
                    rows = [f"{c} ({i})" for i, c in enumerate(["-"]+list(s2))]
                    st.dataframe(pd.DataFrame(mat, columns=cols, index=rows), use_container_width=True)
            
            st.divider()
            c1, c2, c3 = st.columns(3)
            c1.metric("Alignment Score", int(mat[-1, -1]))
            c2.metric("Seq 1 Length", len(s1))
            c3.metric("Seq 2 Length", len(s2))
            st.text("Traceback Path:")
            st.code(f"{res1}\n{res2}", language="text")