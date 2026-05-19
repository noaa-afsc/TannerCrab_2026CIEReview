echo on
copy ..\gmacs.exe                          gmacs.exe
cp ..\gmacs_26_22_03d5_aEffExp1.par        gmacs_26_22_03d5_aEffExp1.par
cp ..\gmacs_26_22_03d5_aEffExp1.Run1.dat   gmacs.dat
cp ..\TannerCrab_Data202509.EffExp1.dat    TannerCrab_Data202509.EffExp1.dat
cp ..\TannerCrab_26_22_03d5_a.Run1.ctl     TannerCrab_26_22_03d5_a.Run1.ctl
cp ..\TannerCrab_26_22_03d5_a.prj          TannerCrab_26_22_03d5_a.prj
gmacs  -rs -phase 8 -nohess -iprint 1 -display 2 -pin gmacs_26_22_03d5_aEffExp1.par -verbose 0 -stopAfterFnCall 1 $@ > chk.rep
