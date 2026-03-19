# 🌿 SPARK Pipeline
### Natural Products Analysis — Visual Interface

---

## What is this?

SPARK is a tool for analyzing natural compounds from plants and other organisms.
You choose a plant family or genus, and the pipeline automatically:
- Extracts all known compounds from the LOTUS database
- Generates statistical charts and chemical diversity figures
- Identifies which compounds are rare and potentially bioactive
- Exports Excel tables and publication-quality PDF figures

---

## Before you start — required software

You need to install **4 programs** before using SPARK.
All of them are free. Follow each step carefully.

---

### Step 1 — Install Python

1. Go to: **https://www.python.org/downloads/**
2. Click the big yellow **"Download Python 3.x.x"** button
3. Run the installer
4. ⚠️ **IMPORTANT:** On the first screen, check the box **"Add Python to PATH"**
   (it is unchecked by default — you must check it)
5. Click **"Install Now"**

![Add Python to PATH checkbox must be checked]

---

### Step 2 — Install R

1. Go to: **https://cran.r-project.org/bin/windows/base/**
2. Click **"Download R x.x.x for Windows"**
3. Run the installer with default options (just keep clicking Next/OK)
4. At the end, click Finish

---

### Step 3 — Install MongoDB

MongoDB is the database that stores the LOTUS compound data.

1. Go to: **https://www.mongodb.com/try/download/community**
2. Under **"MongoDB Community Server"**, select:
   - Version: latest
   - Platform: Windows
   - Package: **msi**
3. Click **Download** and run the installer
4. Keep all default options
5. ✅ Make sure **"Install MongoDB as a Service"** is checked (it is by default)
6. Click through to finish

---

### Step 4 — Install MongoDB Database Tools

This is a **separate package** from MongoDB — it is required to import the LOTUS database.

1. Go to: **https://www.mongodb.com/try/download/database-tools**
2. Select:
   - Platform: Windows
   - Package: **msi**
3. Click **Download** and run the installer with default options

---

### Step 5 (Optional) — Install Rtools

Only needed if you want to use **Part II** (2D structure images and PDF catalog).
If you skip this, everything else still works.

1. Go to: **https://cran.r-project.org/bin/windows/Rtools/**
2. Download the version matching your R version (shown on that page)
3. Run the installer with default options

---

## Getting SPARK

1. Download the ZIP file from the project page
2. Right-click the ZIP → **"Extract All..."**
3. Choose a folder (e.g. `C:\Users\YourName\Documents\SPARK`)
4. Click **Extract**

> ⚠️ Do not extract to a path with special characters (accents, `ç`, spaces in folder names may cause issues on some systems). A simple path like `C:\SPARK` works best.

---

## Getting the database files

You need two types of files placed inside the `DBs\` folder:

### LOTUS database (required for Part I)

Place these files inside `DBs\`:

```
DBs\
├── lotusUniqueNaturalProduct.bson        ← main database (~3.5 GB)
├── lotusUniqueNaturalProduct.metadata.json
├── fragment.bson
├── fragment.metadata.json
├── pubFingerprintsCounts.bson
└── pubFingerprintsCounts.metadata.json
```

### WFO Taxonomic file (recommended)

Place this file inside `DBs\`:

```
DBs\
└── classification.tsv                    ← taxonomic backbone (~900 MB)
```

This file enables taxonomic name normalization (synonyms merged automatically).
Without it, the pipeline still runs but names may be inconsistent.

---

## Running SPARK for the first time

### 1. Double-click `launch.bat`

A black console window will open. It will automatically:

```
[1/6] Checking Python...          ← verifies Python is installed
[2/6] Setting up virtual environment...  ← creates an isolated Python environment
[3/6] Installing Python dependencies...  ← installs all required Python packages
[4/6] Checking R...               ← verifies R is installed
      Checking R packages...      ← installs missing R packages automatically
[5/6] Checking MongoDB...         ← starts the MongoDB database service
[6/6] Launching interface...      ← opens the browser automatically
```

> The **first run takes longer** (5–15 min) because it installs R and Python packages.
> Subsequent runs are fast (under 30 seconds).

> **Keep the black console window open** while using SPARK.
> Closing it stops the application.

---

### 2. The browser opens automatically

Your default browser will open at **http://localhost:8501**

If it doesn't open, manually navigate to that address.

---

### 3. Import the LOTUS database (first time only)

Before running the pipeline, you need to load the LOTUS database into MongoDB.
**This only needs to be done once.**

1. Click the **🗄️ Setup** tab
2. Scroll down to **"Import LOTUS Database into MongoDB"**
3. Confirm the `DBs\` folder path shown is correct
4. Check that the `.bson` files are shown as ✅ found
5. Click **📥 Import into MongoDB**
6. Wait — this takes **15–30 minutes** (3.5 GB of data)
7. When done, you will see ✅ Import finished

> After importing once, you never need to do this again.
> The data stays in MongoDB even after you close SPARK.

---

### 4. Install R packages (first time only)

1. Still in the **🗄️ Setup** tab
2. Scroll up to **"R Packages"**
3. Click **🔍 Check installed packages**
4. If any are missing, click **📥 Install missing packages**
5. Wait for installation to finish (may take 5–15 minutes)

> `launch.bat` already tries to install R packages automatically,
> but you can use this button to verify or retry if something failed.

---

## Running the pipeline

### 5. Configure your analysis

In the left sidebar:

| Setting | What it does |
|---|---|
| **Search level** | Whether you're searching by Family, Genus, or Species |
| **Target taxon(s)** | The name of the family/genus you want to analyze (e.g. `Phyllanthaceae`) |
| **Analysis level** | What to compare in the charts (e.g. genus-level comparison within a family) |

### 6. Choose which modules to run

| Module | What it does | Requires |
|---|---|---|
| **Part I — MongoDB Extraction** | Extracts all compounds for your taxon | MongoDB + LOTUS imported |
| **Part II — Structure Catalog** | Generates 2D structure images and PDF | ChemmineR (Rtools) |
| **Part III — Statistics & Figures** | Statistical plots, heatmaps, diversity charts | Part I cache |
| **Part IV — Bioactivity & Context** | ChEMBL bioactivity + global rarity scores | Part I cache + internet |

### 7. Click ▶️ Start Pipeline

The log window shows real-time progress.
When finished, a green ✅ message appears with the results folder path.

---

## Viewing results

### 📊 Preview tab
- **Metrics:** total compounds, unique structures, taxa count
- **Charts:** interactive bar charts and histograms (hover for values, zoom, download)
- **Image Gallery:** all 2D structure images in a grid
- **Data Tables:** browse any output table directly in the browser

### 📁 Files tab
Download any individual output file.

---

## Output files

Results are saved in `results\lotus_<taxon>_<date>\`:

| File | Contents |
|---|---|
| `*_lin_enriched.parquet` | Full extraction data (used for re-runs without MongoDB) |
| `*_uni_enriched.parquet` | Unique compounds table |
| `*_MASTER_LIST_GLOBAL.xlsx` | Ranked compound list with priority scores |
| `*_BIO_A_Summary.xlsx` | ChEMBL bioactivity data |
| `*_BIO_C_Global_Context.xlsx` | Global rarity per compound |
| `*_Matrix_NatureStyle_V4.pdf` | Novelty × Bioactivity prioritization matrix |
| `*_Profiler_V16_Hydroalcoholic.pdf` | Extraction solvent profiler |
| `png\` | 2D structure images (one per compound) |

---

## Re-running the same taxon

After the first run, **Part I does not need to run again** for the same taxon/date.
The `.parquet` files serve as a cache — Parts II, III and IV load from them instantly.

Simply **uncheck Part I** before clicking Start Pipeline on subsequent runs.

---

## Stopping SPARK

Close the black console window, or press **Ctrl+C** inside it.

To restart, double-click `launch.bat` again.

---

## Common problems

### "Python not found"
You forgot to check **"Add Python to PATH"** during installation.
Uninstall Python and reinstall it, making sure to check that box.

### "R not found"
Add R to the system PATH manually:
1. Press `Win + R`, type `sysdm.cpl`, press Enter
2. Advanced → Environment Variables
3. Under System Variables, find `Path` → Edit → New
4. Add: `C:\Program Files\R\R-4.x.x\bin`
5. Click OK, restart `launch.bat`

### "MongoDB not found" / Parts I and IV-C fail
Make sure MongoDB Community Server is installed (Step 3 above).
The service should start automatically. If it doesn't:
1. Press `Win + R`, type `services.msc`
2. Find **MongoDB** → right-click → Start

### "mongorestore not found" — cannot import database
Install MongoDB Database Tools (Step 4 above).

### "Parquet not found" error when running Parts II/III/IV
Run **Part I** first with MongoDB running and the LOTUS database imported.
This generates the `.parquet` cache files that the other parts need.

### Part II skipped — ChemmineR not installed
Install Rtools (Step 5 above), then use the Setup tab to install missing packages.
Everything else (Parts I, III, IV) works without ChemmineR.

### The window flashes and closes immediately
This usually means Python is not on PATH.
Reinstall Python and check **"Add Python to PATH"**.

### Black console window closes after finishing
This is normal — SPARK stopped.
Double-click `launch.bat` to start again.

---

## Folder structure

```
SPARK\
├── launch.bat                 ← double-click this to start everything
│
├── R\                         ← R analysis scripts (do not modify)
├── app\                       ← interface code (do not modify)
├── DBs\                       ← put your database files here
├── results\                   ← output files appear here after running
└── .venv\                     ← Python environment (created automatically)
```
