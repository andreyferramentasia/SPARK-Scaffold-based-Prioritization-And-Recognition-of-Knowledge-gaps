@echo off
title SPARK-Scaffold-based-Prioritization-And-Recognition-of-Knowledge-gaps
setlocal enabledelayedexpansion

echo.
echo  ============================================================
echo   SPARK-Scaffold-based-Prioritization-And-Recognition-of-Knowledge-gaps - Starting...
echo  ============================================================
echo.

pushd "%~dp0"

set "VENV_PY=.venv\Scripts\python.exe"
set "APP_MAIN=app\main.py"
set "REQUIREMENTS=app\requirements.txt"
set "R_INSTALL=R\install_packages.R"

REM -- 1. Check Python
echo [1/6] Checking Python...
python --version >nul 2>&1
if errorlevel 1 (
    echo.
    echo  ERROR: Python not found.
    echo  Install Python 3.9+ from: https://www.python.org/downloads/
    echo  Check "Add Python to PATH" during installation.
    echo.
    popd
    pause
    exit /b 1
)
for /f "tokens=2" %%v in ('python --version 2^>^&1') do set PY_VER=%%v
echo  OK - Python %PY_VER%

REM -- 2. Create virtual environment
echo.
echo [2/6] Setting up virtual environment...
if not exist "%VENV_PY%" (
    echo  Creating .venv ^(first time only^)...
    python -m venv .venv
    if errorlevel 1 (
        echo  ERROR creating virtual environment.
        popd
        pause
        exit /b 1
    )
    echo  OK - .venv created
) else (
    echo  OK - .venv already exists
)

REM -- 3. Install Python dependencies
echo.
echo [3/6] Installing Python dependencies...
"%VENV_PY%" -m pip install --quiet --upgrade pip >nul 2>&1
"%VENV_PY%" -m pip install --quiet --upgrade -r "%REQUIREMENTS%"
if errorlevel 1 (
    echo  ERROR installing Python packages. Check your internet connection.
    popd
    pause
    exit /b 1
)
echo  OK - Python packages ready

REM -- 4. Check R
echo.
echo [4/6] Checking R...
set "RSCRIPT_EXE="
where Rscript >nul 2>&1
if not errorlevel 1 (
    set "RSCRIPT_EXE=Rscript"
) else (
    for /d %%v in ("C:\Program Files\R\R-*") do (
        if exist "%%v\bin\Rscript.exe" set "RSCRIPT_EXE=%%v\bin\Rscript.exe"
    )
)

if not defined RSCRIPT_EXE (
    echo  WARNING: R not found. Install from: https://www.r-project.org/
    echo  Skipping R package check - use Setup tab after installing R.
    goto skip_r
)

echo  OK - R found: %RSCRIPT_EXE%
echo.
echo  Checking R packages...
"%RSCRIPT_EXE%" --vanilla -e "pkgs=c('dplyr','tidyr','stringr','stringi','data.table','readr','readxl','tibble','forcats','rlang','writexl','jsonlite','arrow','mongolite','httr','progress','ggplot2','ggrepel','scales','viridisLite','vegan','png','ragg','grid','circlize');miss=setdiff(pkgs,rownames(installed.packages()));if(length(miss)){cat('MISSING:',paste(miss,collapse=','),'\n');quit(status=1)}" 2>nul
if errorlevel 1 (
    echo  Installing missing CRAN packages - this may take several minutes...
    "%RSCRIPT_EXE%" --vanilla "%R_INSTALL%"
    echo  Done. Use Setup tab to verify.
) else (
    echo  OK - all CRAN packages installed
)
echo  Bioconductor (ComplexHeatmap, ChemmineR) - optional for Part II...
"%RSCRIPT_EXE%" --vanilla -e "miss=setdiff(c('ComplexHeatmap','ChemmineR'),rownames(installed.packages()));if(length(miss))cat('  Not installed (optional):',paste(miss,collapse=','),'\n') else cat('  OK\n')" 2>nul

:skip_r

REM -- 5. MongoDB
echo.
echo [5/6] Checking MongoDB...
set "MONGOD_EXE="
where mongod >nul 2>&1
if not errorlevel 1 (
    set "MONGOD_EXE=mongod"
) else (
    for /d %%v in ("C:\Program Files\MongoDB\Server\*") do (
        if exist "%%v\bin\mongod.exe" set "MONGOD_EXE=%%v\bin\mongod.exe"
    )
)
if defined MONGOD_EXE (
    netstat -ano 2>nul | find "27017" | find "LISTENING" >nul 2>&1
    if errorlevel 1 (
        echo  Starting MongoDB...
        if not exist "data\db" mkdir "data\db"
        start /B "%MONGOD_EXE%" --dbpath "data\db" --port 27017 --logpath "data\mongod.log" --logappend
        timeout /t 3 /nobreak >nul
        echo  OK - MongoDB started
    ) else (
        echo  OK - MongoDB already running
    )
) else (
    echo  WARNING: MongoDB not found - Parts I and IV-C require it.
    echo  Install from: https://www.mongodb.com/try/download/community
)

REM -- 6. Launch Streamlit
echo.
echo [6/6] Launching interface...
echo.
echo  ============================================================
echo   Open: http://localhost:8501
echo   TO STOP: close this window or press Ctrl+C
echo  ============================================================
echo.

for /f "tokens=5" %%p in ('netstat -ano 2^>nul ^| find "8501" ^| find "LISTENING"') do (
    taskkill /PID %%p /F >nul 2>&1
)
timeout /t 1 /nobreak >nul

start /B cmd /c "timeout /t 4 /nobreak >nul && start http://localhost:8501"

"%VENV_PY%" -m streamlit run "%APP_MAIN%" ^
    --server.headless true ^
    --server.port 8501 ^
    --browser.gatherUsageStats false ^
    --theme.base dark

echo.
echo  Interface closed.
popd
pause
