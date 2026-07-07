@echo off
echo =======================================================================
#          Paper 2 Reproducibility Pipeline Automated Windows Launcher
echo =======================================================================
echo.
echo Checking for Conda installation...
where conda >nul 2>nul
if %errorlevel% neq 0 (
    echo [ERROR] Conda was not found in your PATH. 
    echo Please install Miniconda/Anaconda and make sure to run this script 
    echo from the "Anaconda Prompt".
    pause
    exit /b 1
)

echo Conda found. Creating environment 'amr_env' (this may take a few minutes)...
call conda env create -f environment/environment.yml -y

echo Activating environment 'amr_env'...
call conda activate amr_env

echo Running Paper 2 Reproducibility Pipeline...
python scripts/06_generate_paper2_plots.py

echo.
echo Pipeline execution completed.
pause
