"""
SPARK-Scaffold-based-Prioritization-And-Recognition-of-Knowledge-gaps — Streamlit Interface
Runs R modules (Part I-IV) with a visual configuration panel.
"""

import glob as _glob
import json
import shutil
import socket
import subprocess
import time
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

# ── Paths ─────────────────────────────────────────────────────────────────────
APP_DIR = Path(__file__).parent
PROJECT_DIR = APP_DIR.parent
R_DIR = PROJECT_DIR / "R"
DBS_DIR = PROJECT_DIR / "DBs"
RESULTS_DIR = PROJECT_DIR / "results"
CONFIG_FILE = APP_DIR / "pipeline_config.json"

MAIN_PIPELINE = R_DIR / "Main_Pipeline_new_v2.R"
WFO_TSV = DBS_DIR / "classification.tsv"
WFO_CSV = DBS_DIR / "classification.csv"

# ── Helpers ───────────────────────────────────────────────────────────────────


def find_mongorestore() -> str | None:
    r = shutil.which("mongorestore")
    if r:
        return r
    for pat in [
        r"C:\Program Files\MongoDB\Tools\*\bin\mongorestore.exe",
        r"C:\Program Files\MongoDB\Server\*\bin\mongorestore.exe",
        r"C:\mongodb\bin\mongorestore.exe",
    ]:
        hits = sorted(_glob.glob(pat))
        if hits:
            return hits[-1]
    return None


def import_bson(mongorestore: str, bson_path: Path, db: str, collection: str) -> tuple[bool, str]:
    cmd = [
        mongorestore,
        "--uri", "mongodb://127.0.0.1:27017",
        "--db", db,
        "--collection", collection,
        "--drop",
        str(bson_path),
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           encoding="utf-8", errors="replace", timeout=3600)
        if r.returncode == 0:
            return True, "OK"
        return False, r.stderr[-500:] or r.stdout[-500:]
    except subprocess.TimeoutExpired:
        return False, "Timeout after 1 hour."
    except Exception as e:
        return False, str(e)


def find_rscript() -> str | None:
    r = shutil.which("Rscript")
    if r:
        return r
    for pat in [
        r"C:\Program Files\R\R-*\bin\Rscript.exe",
        r"C:\Program Files (x86)\R\R-*\bin\Rscript.exe",
    ]:
        hits = sorted(_glob.glob(pat))
        if hits:
            return hits[-1]
    return None


def is_mongo_running() -> bool:
    try:
        with socket.create_connection(("127.0.0.1", 27017), timeout=1):
            return True
    except OSError:
        return False


def wfo_path() -> str | None:
    if WFO_TSV.exists():
        return str(WFO_TSV)
    if WFO_CSV.exists():
        return str(WFO_CSV)
    return None


def latest_result_dir(out_dir: Path, tag_prefix: str) -> Path | None:
    if not out_dir.exists():
        return None
    matches = sorted(
        [d for d in out_dir.iterdir() if d.is_dir(
        ) and d.name.startswith(f"lotus_{tag_prefix}")],
        key=lambda d: d.stat().st_mtime,
        reverse=True,
    )
    return matches[0] if matches else None


def list_result_files(folder: Path) -> list[Path]:
    return sorted(
        [f for f in folder.rglob("*") if f.is_file()],
        key=lambda f: f.stat().st_size,
        reverse=True,
    )[:40]


def build_tag_prefix(mode: str, values: list[str]) -> str:
    def safe(s): return "".join(c if c.isalnum()
                                or c in "._-" else "_" for c in s)
    return safe(f"{mode}_{'_'.join(values)}")


# ── Page setup ────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="SPARK-Scaffold-based-Prioritization-And-Recognition-of-Knowledge-gaps",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    .block-container { padding-top: 1.5rem; }
    div[data-testid="stVerticalBlock"] > div { gap: 0.4rem; }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("🌿 SPARK-Scaffold-based-Prioritization-And-Recognition-of-Knowledge-gaps")
st.caption("Visual interface for natural products analysis — extraction, structures, statistics and bioactivity.")

# ── Session state defaults ────────────────────────────────────────────────────
if "running" not in st.session_state:
    st.session_state.running = False
if "log_lines" not in st.session_state:
    st.session_state.log_lines = []
if "result_dir" not in st.session_state:
    st.session_state.result_dir = None

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("⚙️ Configuration")

    st.subheader("Taxonomic Target")
    taxon_mode = st.selectbox(
        "Search level",
        ["family", "genus", "species"],
        index=0,
        help="Taxonomic level used to filter the LOTUS database.",
        disabled=st.session_state.running,
    )
    taxon_values_raw = st.text_input(
        "Target taxon(s)",
        value="Phyllanthaceae",
        help="One or more names separated by comma. E.g.: Phyllanthaceae or Ocotea, Nectandra",
        disabled=st.session_state.running,
    )
    taxon_values = [v.strip()
                    for v in taxon_values_raw.split(",") if v.strip()]

    analysis_level = st.selectbox(
        "Analysis level (comparison in plots)",
        ["genus", "species", "family"],
        index=0,
        help="What each column/row represents in the figures.",
        disabled=st.session_state.running,
    )

    st.divider()

    st.subheader("Modules")
    run_p1 = st.checkbox(
        "Part I — MongoDB Extraction",
        value=True,
        help="Extracts compounds from the LOTUS database. Requires MongoDB running.",
        disabled=st.session_state.running,
    )
    run_p2 = st.checkbox(
        "Part II — Structure Catalog (PDF)",
        value=True,
        help="Generates PNG images and PDF catalog of 2D structures. Requires ChemmineR.",
        disabled=st.session_state.running,
    )
    run_p3 = st.checkbox(
        "Part III — Statistics & Figures",
        value=True,
        help="Produces statistical plots and summary tables.",
        disabled=st.session_state.running,
    )
    run_p4 = st.checkbox(
        "Part IV — Bioactivity & Global Context",
        value=True,
        help="Integrates ChEMBL bioactivity data and global rarity context (MongoDB).",
        disabled=st.session_state.running,
    )

    if run_p1 or run_p4:
        st.info("⚠️ Parts I and IV-C require MongoDB on port 27017.")

    st.divider()

    with st.expander("Advanced Settings"):
        top_taxa = st.slider("Max taxa in plots", 5, 100,
                             40, disabled=st.session_state.running)
        min_compounds = st.slider(
            "Min compounds per taxon", 1, 50, 10, disabled=st.session_state.running)
        use_wfo = st.checkbox("WFO normalization", value=True,
                              disabled=st.session_state.running)
        export_excel = st.checkbox(
            "Export Excel (.xlsx)", value=True, disabled=st.session_state.running)
        out_dir_str = st.text_input(
            "Results directory",
            value=str(RESULTS_DIR),
            disabled=st.session_state.running,
        )

    st.divider()

    st.subheader("R Installation")
    detected_r = find_rscript()
    r_path = st.text_input(
        "Rscript executable",
        value=detected_r or "",
        placeholder=r"C:\Program Files\R\R-4.x.x\bin\Rscript.exe",
        help="Full path to Rscript.exe",
        disabled=st.session_state.running,
    )

# ── Status bar ────────────────────────────────────────────────────────────────
col_r, col_mongo, col_wfo = st.columns(3)

r_ok = bool(r_path) and Path(r_path).exists()
with col_r:
    if r_ok:
        st.success("✅ R found")
    else:
        st.error("❌ R not found — check the path")

needs_mongo = run_p1 or run_p4
with col_mongo:
    if needs_mongo:
        mongo_ok = is_mongo_running()
        if mongo_ok:
            st.success("✅ MongoDB connected (port 27017)")
        else:
            st.warning("⚠️ MongoDB not found — Parts I/IV-C may fail")
    else:
        st.info("ℹ️ MongoDB not required for selected modules")

wfo_file = wfo_path()
with col_wfo:
    if wfo_file:
        st.success(f"✅ WFO found: {Path(wfo_file).name}")
    elif use_wfo:
        st.warning("⚠️ classification.tsv not found in DBs/")
    else:
        st.info("ℹ️ WFO disabled")

st.divider()

# ── Tabs ──────────────────────────────────────────────────────────────────────
tab_setup, tab_run, tab_preview, tab_results, tab_help = st.tabs(
    ["🗄️ Setup", "▶️ Run", "📊 Preview", "📁 Files", "❓ Help"]
)

# ── TAB: RUN ──────────────────────────────────────────────────────────────────
with tab_run:
    if not taxon_values:
        st.warning("Enter at least one target taxon in the sidebar.")

    can_run = r_ok and bool(taxon_values) and (
        run_p1 or run_p2 or run_p3 or run_p4)

    run_btn = st.button(
        "▶️ Start Pipeline",
        disabled=not can_run or st.session_state.running,
        type="primary",
        use_container_width=True,
    )

    stop_btn = st.button(
        "⏹ Stop",
        disabled=not st.session_state.running,
        use_container_width=True,
    )

    log_area = st.empty()

    def display_log(lines: list[str]) -> None:
        text = "\n".join(lines[-80:])
        log_area.code(text or "(waiting for R output...)", language=None)

    display_log(st.session_state.log_lines)

    if stop_btn and "proc" in st.session_state:
        try:
            st.session_state.proc.terminate()
        except Exception:
            pass
        st.session_state.running = False
        st.warning("Pipeline stopped by user.")

    if run_btn:
        st.session_state.running = True
        st.session_state.log_lines = []
        st.session_state.result_dir = None

        out_dir = Path(out_dir_str)
        out_dir.mkdir(parents=True, exist_ok=True)

        config = {
            "taxon_mode":                       taxon_mode,
            "taxon_values":                     taxon_values,
            "analysis_tax_level":               analysis_level,
            "run_module1":                      run_p1,
            "run_module2":                      run_p2,
            "run_module3":                      run_p3,
            "run_module4":                      run_p4,
            "analysis_top_taxa":                top_taxa,
            "analysis_min_compounds_per_taxon": min_compounds,
            "out_dir_base":                     str(out_dir),
            "use_WFO_normalization":            use_wfo and bool(wfo_file),
            "wfo_csv_path":                     wfo_file or "",
            "export_excel":                     export_excel,
            "export_parquet":                   True,
            "verbose":                          True,
            "mongo_url":                        "mongodb://127.0.0.1:27017/?socketTimeoutMS=3600000&connectTimeoutMS=300000&serverSelectionTimeoutMS=300000",
            "db_name":                          "lotus",
            "coll_name":                        "lotusUniqueNaturalProduct",
            "r_scripts_dir":                    str(R_DIR),
        }

        CONFIG_FILE.write_text(json.dumps(
            config, indent=2, ensure_ascii=True), encoding="ascii")

        proc = subprocess.Popen(
            [r_path, str(MAIN_PIPELINE), str(CONFIG_FILE)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            cwd=str(R_DIR),
        )
        st.session_state.proc = proc

        ts_start = time.time()
        with st.spinner("Pipeline running — follow the log below..."):
            for line in proc.stdout:
                st.session_state.log_lines.append(line.rstrip())
                display_log(st.session_state.log_lines)
            returncode = proc.wait()

        elapsed = time.time() - ts_start
        st.session_state.running = False

        if returncode == 0:
            st.success(f"✅ Pipeline completed in {elapsed:.0f}s")
            tag_prefix = build_tag_prefix(taxon_mode, taxon_values)
            result = latest_result_dir(out_dir, tag_prefix)
            if result:
                st.session_state.result_dir = result
                st.info(f"📁 Results saved to: `{result}`")
        else:
            st.error(f"❌ Pipeline exited with error code {returncode}")
            st.info("Check the log above for details.")

# ── TAB: SETUP ────────────────────────────────────────────────────────────────
with tab_setup:

    # ── R Packages ─────────────────────────────────────────────────────────────
    st.subheader("📦 R Packages")

    CHECK_SCRIPT = R_DIR / "check_packages.R"

    def check_r_packages(rscript: str) -> dict | None:
        try:
            r = subprocess.run(
                [rscript, "--vanilla", str(CHECK_SCRIPT)],
                capture_output=True, text=True, timeout=30,
                encoding="utf-8", errors="replace"
            )
            return json.loads(r.stdout.strip())
        except Exception:
            return None

    if r_ok:
        if st.button("🔍 Check installed packages", use_container_width=True):
            with st.spinner("Checking..."):
                pkg_status = check_r_packages(r_path)

            if pkg_status is None:
                st.error("Could not check packages. Verify the R path.")
            else:
                missing = [p for p, ok in pkg_status.items() if not ok]
                installed_count = sum(pkg_status.values())
                total = len(pkg_status)
                st.write(f"**{installed_count}/{total} packages installed**")

                cols = st.columns(3)
                for i, (pkg, ok) in enumerate(sorted(pkg_status.items())):
                    cols[i % 3].write(f"{'✅' if ok else '❌'} {pkg}")

                if missing:
                    st.warning(f"Missing: {', '.join(missing)}")
                    if st.button("📥 Install missing packages", type="primary", use_container_width=True):
                        install_script = R_DIR / "install_packages.R"
                        log_ph = st.empty()
                        lines = []
                        proc = subprocess.Popen(
                            [r_path, "--vanilla", str(install_script)],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, encoding="utf-8", errors="replace",
                            cwd=str(R_DIR)
                        )
                        for line in proc.stdout:
                            lines.append(line.rstrip())
                            log_ph.code("\n".join(lines[-40:]))
                        proc.wait()
                        if proc.returncode == 0:
                            st.success(
                                "✅ Installation complete! Check packages again to confirm.")
                        else:
                            st.error(
                                "Installation finished with errors — see the log above.")
                else:
                    st.success("✅ All packages are installed!")
    else:
        st.info("Set the R path in the sidebar to check packages.")

    st.divider()

    # ── MongoDB Import ──────────────────────────────────────────────────────────
    st.subheader("🗄️ Import LOTUS Database into MongoDB")

    mongo_live = is_mongo_running()
    if mongo_live:
        st.success("✅ MongoDB is running on port 27017")
    else:
        st.error(
            "❌ MongoDB is not running. Install and start the service before importing.")
        st.markdown("Download: https://www.mongodb.com/try/download/community")

    st.divider()

    mr = find_mongorestore()
    if mr:
        st.success(f"✅ mongorestore found: `{mr}`")
    else:
        st.error("❌ mongorestore not found.")
        st.markdown(
            "**MongoDB Database Tools** is a separate package. Download at:\n\n"
            "https://www.mongodb.com/try/download/database-tools\n\n"
            "Select **Windows / msi** and install with default options."
        )

    st.divider()

    st.markdown("**BSON files folder:**")
    dbs_input = st.text_input(
        "DBs folder path",
        value=str(DBS_DIR),
        help="Folder containing lotusUniqueNaturalProduct.bson, fragment.bson, etc.",
    )
    dbs_path = Path(dbs_input) if dbs_input else DBS_DIR

    bsons = {
        "lotusUniqueNaturalProduct": dbs_path / "lotusUniqueNaturalProduct.bson",
        "fragment":                  dbs_path / "fragment.bson",
        "pubFingerprintsCounts":     dbs_path / "pubFingerprintsCounts.bson",
    }

    st.markdown("**Files found:**")
    all_found = True
    for col_name, bson_file in bsons.items():
        if bson_file.exists():
            size_mb = bson_file.stat().st_size / 1024 / 1024
            st.write(f"✅ `{bson_file.name}` — {size_mb:,.0f} MB")
        else:
            st.write(f"❌ `{bson_file.name}` — not found")
            if col_name == "lotusUniqueNaturalProduct":
                all_found = False

    st.divider()

    can_import = mongo_live and mr and all_found
    if st.button("📥 Import into MongoDB", disabled=not can_import, type="primary", use_container_width=True):
        st.warning(
            "⏳ Importing... Do not close this window. May take 15–30 min for 3.5 GB.")
        progress = st.progress(0, text="Starting import...")
        log_box = st.empty()
        import_log = []

        total = sum(1 for b in bsons.values() if b.exists())
        done = 0

        for col_name, bson_file in bsons.items():
            if not bson_file.exists():
                continue
            progress.progress(done / total, text=f"Importing {col_name}...")
            import_log.append(
                f"-> {col_name} ({bson_file.stat().st_size/1024/1024:,.0f} MB)...")
            log_box.code("\n".join(import_log))

            ok, msg = import_bson(mr, bson_file, "lotus", col_name)
            import_log.append(f"  {'OK' if ok else 'ERROR: ' + msg}")
            log_box.code("\n".join(import_log))
            done += 1
            progress.progress(done / total, text=f"Done: {col_name}")

        progress.progress(1.0, text="Import finished!")
        st.success("✅ LOTUS database imported! Now run Part I in the Run tab.")

    if not mongo_live:
        st.info("Start MongoDB first to enable import.")
    elif not mr:
        st.info("Install MongoDB Database Tools to enable import.")
    elif not all_found:
        st.info("lotusUniqueNaturalProduct.bson not found in the specified folder.")

# ── TAB: PREVIEW ──────────────────────────────────────────────────────────────
with tab_preview:

    def _find_taxon_col(df: pd.DataFrame, hints: list[str]) -> str | None:
        """Return the first column whose name contains any of the hint keywords."""
        for hint in hints:
            for col in df.columns:
                if hint in col.lower():
                    return col
        return None

    def _load_df(path: Path, nrows: int | None = None) -> pd.DataFrame | None:
        try:
            if path.suffix == ".parquet":
                df = pd.read_parquet(path)
                return df.head(nrows) if nrows else df
            if path.suffix == ".xlsx":
                return pd.read_excel(path, nrows=nrows)
            if path.suffix == ".csv":
                return pd.read_csv(path, nrows=nrows)
        except Exception:
            return None

    out_dir_prev = Path(out_dir_str) if "out_dir_str" in dir() else RESULTS_DIR
    prev_dirs = []
    if out_dir_prev.exists():
        prev_dirs = sorted(
            [d for d in out_dir_prev.iterdir() if d.is_dir()],
            key=lambda d: d.stat().st_mtime,
            reverse=True,
        )

    if not prev_dirs:
        st.info("No results found yet. Run the pipeline first.")
    else:
        sel_prev = st.selectbox(
            "Select a run",
            options=prev_dirs,
            format_func=lambda d: d.name,
            key="preview_sel",
        )
        if sel_prev:
            sel_prev = Path(sel_prev)

            # ── Locate main parquet files ──────────────────────────────────────
            all_pq = list(sel_prev.glob("*.parquet"))
            lin_pq = next(
                (f for f in all_pq if "lin_enriched" in f.name), None)
            uni_pq = next(
                (f for f in all_pq if "uni_enriched" in f.name), None)
            main_df = _load_df(uni_pq) if uni_pq else (
                _load_df(lin_pq) if lin_pq else None)

            # ── Metrics row ───────────────────────────────────────────────────
            if main_df is not None:
                tax_col = _find_taxon_col(
                    main_df, ["genus", "species", "family", "taxon"])
                inchi_col = _find_taxon_col(
                    main_df, ["inchikey", "inchi_key", "structure"])

                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Rows (compounds)", f"{len(main_df):,}")
                m1.caption(Path(uni_pq or lin_pq).name)
                if inchi_col:
                    m2.metric("Unique structures",
                              f"{main_df[inchi_col].nunique():,}")
                if tax_col:
                    m3.metric(
                        f"Unique {tax_col.split('_')[-1]}", f"{main_df[tax_col].nunique():,}")
                total_mb = sum(f.stat().st_size for f in sel_prev.rglob(
                    "*") if f.is_file()) / 1e6
                m4.metric("Total output size", f"{total_mb:.1f} MB")

                st.divider()

            # ── Charts ────────────────────────────────────────────────────────
            chart_df = _load_df(lin_pq) if lin_pq else main_df
            if chart_df is not None:
                st.subheader("📈 Charts")
                ctab1, ctab2, ctab3 = st.tabs(
                    ["Top taxa", "Compound distribution", "Heatmap preview"])

                with ctab1:
                    tax_col = _find_taxon_col(
                        chart_df, ["genus", "species", "family"])
                    if tax_col:
                        top_n = st.slider("Show top N taxa",
                                          5, 50, 20, key="chart_topn")
                        counts = (
                            chart_df[tax_col]
                            .value_counts()
                            .head(top_n)
                            .reset_index()
                        )
                        counts.columns = [tax_col, "count"]
                        fig = px.bar(
                            counts.sort_values("count"),
                            x="count", y=tax_col,
                            orientation="h",
                            title=f"Compounds per {tax_col.split('_')[-1]} (top {top_n})",
                            labels={"count": "# compounds", tax_col: ""},
                            color="count",
                            color_continuous_scale="teal",
                        )
                        fig.update_layout(
                            height=max(300, top_n * 22),
                            showlegend=False,
                            coloraxis_showscale=False,
                            plot_bgcolor="rgba(0,0,0,0)",
                            paper_bgcolor="rgba(0,0,0,0)",
                            font_color="#ccc",
                        )
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        st.info("No taxonomic column detected in the data.")

                with ctab2:
                    num_cols = chart_df.select_dtypes(
                        include="number").columns.tolist()
                    num_cols = [
                        c for c in num_cols if chart_df[c].nunique() > 2]
                    if num_cols:
                        chosen_num = st.selectbox(
                            "Column to plot", num_cols, key="hist_col")
                        fig2 = px.histogram(
                            chart_df,
                            x=chosen_num,
                            nbins=40,
                            title=f"Distribution of {chosen_num}",
                            color_discrete_sequence=["#2ecc71"],
                        )
                        fig2.update_layout(
                            plot_bgcolor="rgba(0,0,0,0)",
                            paper_bgcolor="rgba(0,0,0,0)",
                            font_color="#ccc",
                        )
                        st.plotly_chart(fig2, use_container_width=True)
                    else:
                        st.info("No numeric columns available for histogram.")

                with ctab3:
                    tax_col2 = _find_taxon_col(chart_df, ["genus", "family"])
                    num_col2 = next(
                        (c for c in chart_df.select_dtypes(include="number").columns
                         if chart_df[c].nunique() > 2), None
                    )
                    if tax_col2 and num_col2:
                        top20 = chart_df[tax_col2].value_counts().head(
                            20).index
                        pivot = (
                            chart_df[chart_df[tax_col2].isin(top20)]
                            .groupby(tax_col2)[num_col2]
                            .mean()
                            .reset_index()
                            .sort_values(num_col2, ascending=False)
                        )
                        fig3 = px.bar(
                            pivot, x=tax_col2, y=num_col2,
                            title=f"Mean {num_col2} by {tax_col2.split('_')[-1]} (top 20)",
                            color=num_col2,
                            color_continuous_scale="viridis",
                        )
                        fig3.update_layout(
                            plot_bgcolor="rgba(0,0,0,0)",
                            paper_bgcolor="rgba(0,0,0,0)",
                            font_color="#ccc",
                        )
                        st.plotly_chart(fig3, use_container_width=True)
                    else:
                        st.info(
                            "Need a taxon column + numeric column to build this chart.")

                st.divider()

            # ── Image Gallery ─────────────────────────────────────────────────
            images = sorted(sel_prev.rglob("*.png")) + \
                sorted(sel_prev.rglob("*.jpg"))
            images = [i for i in images if i.stat().st_size > 2048][:80]

            if images:
                st.subheader(f"🖼️ Image Gallery ({len(images)})")

                # Filter by subfolder
                subfolders = sorted({i.parent.name for i in images})
                if len(subfolders) > 1:
                    folder_filter = st.selectbox(
                        "Filter by folder", ["All"] + subfolders, key="img_filter"
                    )
                    if folder_filter != "All":
                        images = [
                            i for i in images if i.parent.name == folder_filter]

                cols_per_row = st.select_slider(
                    "Columns", options=[2, 3, 4, 5], value=3, key="img_cols"
                )
                rows = [images[i:i+cols_per_row]
                        for i in range(0, len(images), cols_per_row)]
                for row in rows:
                    img_cols = st.columns(cols_per_row)
                    for col, img_path in zip(img_cols, row):
                        col.image(str(img_path), caption=img_path.name,
                                  use_container_width=True)

                st.divider()

            # ── Data Tables ───────────────────────────────────────────────────
            data_files = (
                sorted(sel_prev.rglob("*.parquet"))
                + sorted(sel_prev.rglob("*.xlsx"))
                + sorted(sel_prev.rglob("*.csv"))
            )
            data_files = [f for f in data_files if f.stat().st_size >
                          1024][:12]

            if data_files:
                st.subheader(f"📋 Data Tables ({len(data_files)})")
                chosen = st.selectbox(
                    "Select table",
                    options=data_files,
                    format_func=lambda f: f.relative_to(sel_prev).as_posix(),
                    key="preview_table",
                )
                if chosen:
                    df_t = _load_df(Path(chosen), nrows=500)
                    if df_t is not None:
                        st.caption(
                            f"{len(df_t)} rows × {len(df_t.columns)} columns")
                        st.dataframe(
                            df_t, use_container_width=True, height=380)
                    else:
                        st.error("Could not read this file.")

            if not images and chart_df is None and not data_files:
                st.info("No output files found in this run folder.")

# ── TAB: RESULTS ──────────────────────────────────────────────────────────────
with tab_results:
    out_dir_p = Path(out_dir_str) if "out_dir_str" in dir() else RESULTS_DIR
    result_dirs = []
    if out_dir_p.exists():
        result_dirs = sorted(
            [d for d in out_dir_p.iterdir() if d.is_dir()],
            key=lambda d: d.stat().st_mtime,
            reverse=True,
        )

    if not result_dirs:
        st.info("No results found yet. Run the pipeline first.")
    else:
        selected = st.selectbox(
            "Select a run",
            options=result_dirs,
            format_func=lambda d: d.name,
        )
        if selected:
            files = list_result_files(Path(selected))
            if not files:
                st.info("Empty folder.")
            else:
                st.write(f"**{len(files)} file(s) — click to download:**")
                for fpath in files:
                    rel = fpath.relative_to(selected)
                    size_kb = fpath.stat().st_size / 1024
                    with open(fpath, "rb") as fh:
                        st.download_button(
                            label=f"⬇️ {rel}  ({size_kb:,.0f} KB)",
                            data=fh.read(),
                            file_name=fpath.name,
                            key=str(fpath),
                            use_container_width=True,
                        )

# ── TAB: HELP ─────────────────────────────────────────────────────────────────
with tab_help:
    st.markdown(
        "## How to use SPARK-Scaffold-based-Prioritization-And-Recognition-of-Knowledge-gaps")

    st.markdown("""
**Recommended order on first run:**
1. 🗄️ **Setup** — install R packages, import MongoDB database
2. ▶️ **Run** — run Part I (extraction) first, then Parts II–IV
3. 📊 **Preview** — explore charts, images and tables
4. 📁 **Files** — download all output files
""")

    st.markdown("""
```
Part I   → connects to MongoDB → extracts compounds → saves lin_enriched.parquet + uni_enriched.parquet
Part II  → reads parquet → renders 2D structure images → generates PDF catalog  (optional, needs ChemmineR)
Part III → reads parquet → produces statistical plots + summary tables
Part IV  → reads parquet + MongoDB → ChEMBL bioactivity + global rarity context
```
> **Part I only needs to run once per taxon/date.** Parts II–IV reload from the Parquet cache automatically.
""")

    st.divider()
    st.subheader("🔧 Troubleshooting")

    with st.expander("❌  Python not found / 'python' is not recognized"):
        st.markdown("""
**Cause:** Python is not installed or not on the system PATH.

**Fix:**
1. Download Python 3.9+ from https://www.python.org/downloads/
2. During installation, check **"Add Python to PATH"** (critical!)
3. Restart the terminal / launch.bat and try again

**Verify:** open CMD and run `python --version` — should show `Python 3.x.x`.
""")

    with st.expander("❌  R not found / 'Rscript' is not recognized"):
        st.markdown("""
**Cause:** R is not installed or its `bin/` folder is not on PATH.

**Fix:**
1. Download R 4.x from https://www.r-project.org/
2. During installation, allow it to register on PATH (default)
3. Or manually add `C:\\Program Files\\R\\R-4.x.x\\bin` to your system PATH
4. Restart launch.bat

**Verify:** open CMD and run `Rscript --version`.
""")

    with st.expander("❌  R package not found / 'there is no package called ...'"):
        st.markdown("""
**Cause:** One or more R packages are not installed in your R library.

**Fix — via interface:**
1. Go to **Setup** tab → R Packages → click **Check installed packages**
2. If packages are missing, click **Install missing packages**
3. Wait for installation to finish (may take 5–15 min on first run)

**Fix — via R console:**
```r
source("R/install_packages.R")
```

**Common missing packages:** `arrow`, `mongolite`, `ggrepel`, `vegan`, `ragg`, `circlize`
""")

    with st.expander("❌  ChemmineR / ComplexHeatmap not installed (Part II fails)"):
        st.markdown("""
**Cause:** ChemmineR and ComplexHeatmap are Bioconductor packages that require
**Rtools** on Windows to compile from source.

**Part II is optional** — if these packages are missing, Part II is automatically
skipped and Parts I, III, IV continue normally.

**To install ChemmineR (optional):**
1. Download Rtools from https://cran.r-project.org/bin/windows/Rtools/
2. Install Rtools with default options (adds compiler to PATH)
3. Restart your machine
4. In the Setup tab, click **Install missing packages**

Or in R console:
```r
BiocManager::install("ChemmineR")
BiocManager::install("ComplexHeatmap")
```
""")

    with st.expander("❌  MongoDB not running / connection refused on port 27017"):
        st.markdown("""
**Cause:** The MongoDB service is not started.

**Fix — automatic:** launch.bat tries to start MongoDB automatically if it finds
`mongod.exe` in the standard installation path.

**Fix — manual:**
1. Open CMD as Administrator
2. Run: `net start MongoDB`
   or: `"C:\\Program Files\\MongoDB\\Server\\8.0\\bin\\mongod.exe" --dbpath C:\\data\\db`

**MongoDB not installed?**
1. Download from https://www.mongodb.com/try/download/community
2. Install as a **Windows Service** (default option)
3. The service starts automatically on boot after that

**Verify:** open CMD and run `netstat -ano | find "27017"` — should show LISTENING.
""")

    with st.expander("❌  mongorestore not found (cannot import BSON)"):
        st.markdown("""
**Cause:** `mongorestore` is part of **MongoDB Database Tools**, which is a
**separate download** from MongoDB Server.

**Fix:**
1. Download from https://www.mongodb.com/try/download/database-tools
2. Select **Windows / msi** and install with default options
3. The tools are added to PATH automatically
4. Reopen launch.bat

**Verify:** open CMD and run `mongorestore --version`.
""")

    with st.expander("❌  LOTUS database not imported / 'collection is empty'"):
        st.markdown("""
**Cause:** MongoDB is running but the LOTUS data has not been imported yet.

**Fix:**
1. Make sure you have the `.bson` files in the `DBs/` folder:
   - `lotusUniqueNaturalProduct.bson` (required, ~3.5 GB)
   - `fragment.bson` (optional)
   - `pubFingerprintsCounts.bson` (optional)
2. Go to **Setup** tab → Import LOTUS Database
3. Click **Import into MongoDB** and wait (15–30 min for 3.5 GB)

> You only need to import **once**. After that the data persists in MongoDB.

**Verify in MongoDB Compass:** connect to `localhost:27017` → look for database `lotus`
→ collection `lotusUniqueNaturalProduct` — should have millions of documents.
""")

    with st.expander("❌  Parquet not found / 'ERRO CRÍTICO: Arquivos Parquet não encontrados'"):
        st.markdown("""
**Cause:** Parts II–IV require the `.parquet` cache files generated by Part I,
but Part I hasn't run yet (or failed before saving).

**Fix:**
1. Check **Part I — MongoDB Extraction** in the Run tab
2. Make sure MongoDB is running and the LOTUS database is imported
3. Run the pipeline — Part I will generate the cache files
4. After that you can run Parts II–IV without Part I

**Where are the files?**
`results/lotus_<taxon>_<date>/lotus_<taxon>_<date>_lin_enriched.parquet`
`results/lotus_<taxon>_<date>/lotus_<taxon>_<date>_uni_enriched.parquet`
""")

    with st.expander("❌  WFO file not found / taxonomic normalization disabled"):
        st.markdown("""
**Cause:** The WFO (World Flora Online) taxonomic backbone file is missing from `DBs/`.

**Expected file:** `DBs/classification.tsv` (~900 MB)

**Impact:** Taxonomic normalization is skipped — results may have inconsistent
taxon names (synonyms not merged).

**Fix:**
- Download `classification.tsv` from the World Flora Online backbone release
- Place it in the `DBs/` folder
- Or **uncheck WFO normalization** in Advanced Settings to run without it
""")

    with st.expander("❌  Pipeline exits with error code / generic R error in the log"):
        st.markdown("""
**Steps to diagnose:**
1. Read the full log in the **Run** tab — the last error line usually identifies the problem
2. Common patterns:
   - `Error in library(X)` → missing R package → go to Setup → install packages
   - `could not find function` → wrong R version or package version mismatch
   - `Error: MONGO_ERR` → MongoDB not running or data not imported
   - `subscript out of bounds` / `object not found` → data issue, check taxon name spelling
   - `cannot allocate vector of size X` → out of RAM (close other programs and retry)

**Taxon name issues:**
- Names are case-sensitive in the LOTUS database
- Use the exact family/genus name as it appears in the database
- Example: `Phyllanthaceae` not `phyllanthaceae`

**Still stuck?** Copy the last 10 lines from the log and search online or check the LOTUS database documentation.
""")

    with st.expander("❌  Port 8501 already in use"):
        st.markdown("""
**Cause:** A previous Streamlit instance is still running.

**Fix — automatic:** launch.bat kills any process on port 8501 before starting.

**Fix — manual:** open CMD and run:
```
for /f "tokens=5" %p in ('netstat -ano ^| find "8501" ^| find "LISTENING"') do taskkill /PID %p /F
```
Then rerun launch.bat.
""")

    with st.expander("❌  Interface opens but shows a blank page / old content"):
        st.markdown("""
**Cause:** Browser cached an old version of the app.

**Fix:** Press **Ctrl+Shift+R** (hard reload) in the browser.

If the wrong app still shows, check if multiple Streamlit processes are running:
open CMD → `netstat -ano | find "8501"` — should show only one PID.
""")

    st.divider()
    st.subheader("📂 Required files reference")
    st.markdown("""
| File | Location | Required for |
|---|---|---|
| `lotusUniqueNaturalProduct.bson` | `DBs/` | Part I (MongoDB import) |
| `fragment.bson` | `DBs/` | Part I (optional) |
| `pubFingerprintsCounts.bson` | `DBs/` | Part I (optional) |
| `classification.tsv` | `DBs/` | WFO normalization (optional) |
| `Main_Pipeline_new_v2.R` | `R/` | All parts |
| `Part I - Extraction_new_v2.R` | `R/` | Part I |
| `Part II - Structures_new_v2.R` | `R/` | Part II (optional) |
| `Part III - Figure_Statistic_new_v2.R` | `R/` | Part III |
| `Part IV - actives_occurence_v3.R` | `R/` | Part IV |

**Output location:** `results/lotus_<taxon>_<analysis_level>_<date>/`
- `*.parquet` — data cache (re-used by subsequent runs)
- `*.xlsx` — Excel tables
- `*.pdf` — figures and structure catalog
- `png/` — 2D structure images
- `diagnostics/` — QC logs
""")

    st.subheader("⚡ Performance tips")
    st.markdown("""
- **Part I is slow (~20–60 min)** on first run — it queries 3.5 GB from MongoDB. Let it finish.
- After Part I completes, **uncheck Part I** for subsequent runs — data loads from the Parquet cache instantly.
- **Part II** generates one PNG per compound — can take a long time for large datasets.
- If you only need statistics, run **Parts III and IV** with Part I cached.
- MongoDB uses ~1–2 GB RAM while running — close other heavy programs.
""")
