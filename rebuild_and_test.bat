C:
CD C:\users\szelazni\desktop\github\black_scholes

python2 setup.py build_ext
IF %ERRORLEVEL% GEQ 1 GOTO HANDLER

python2 setup.py install
IF %ERRORLEVEL% GEQ 1 GOTO HANDLER

python3 setup.py build_ext
IF %ERRORLEVEL% GEQ 1 GOTO HANDLER

python3 setup.py install
IF %ERRORLEVEL% GEQ 1 GOTO HANDLER

python2 black_scholes__unittest.py -v
python3 black_scholes__unittest.py -v

:HANDLER
pause
