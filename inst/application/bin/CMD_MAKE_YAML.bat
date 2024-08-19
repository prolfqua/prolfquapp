@echo off

REM Get the path to the installed R package
for /f "tokens=*" %%i in ('Rscript --vanilla -e "cat(system.file(package = 'prolfquapp'))"') do set PACKAGE_PATH=%%i

REM Build the path to the R script
set R_SCRIPT_PATH=%PACKAGE_PATH%\application\CMD_MAKE_YAML.R

REM Check if the R script exists
if exist "%R_SCRIPT_PATH%" (
    Rscript --vanilla "%R_SCRIPT_PATH%" %*
) else (
    echo Error: R script not found at '%R_SCRIPT_PATH%'
    exit /b 1
)


