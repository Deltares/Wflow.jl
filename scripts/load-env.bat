@echo off
if exist "%PIXI_PROJECT_ROOT%\.env" (
    for /f "delims=" %%A in ('dotenv -f "%PIXI_PROJECT_ROOT%\.env" list --format=simple') do set "%%A"
)
