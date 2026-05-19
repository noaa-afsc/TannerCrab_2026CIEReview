echo on
copy ..\tcsam02.exe tcsam02.exe
tcsam02  -rs -nox  -configFile ../MCI.2025.22_03d5.inp -phase 5  -calcOFL  -binp ../tcsam02.best.bar 
tcsam02  -rs -nox  -configFile ../MCI.2025.22_03d5.inp -phase 5  -calcOFL  -binp ../tcsam02.bar -hess_step 5


