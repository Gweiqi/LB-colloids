SET PATH=%PATH%;C:\TDM-GCC-64\bin
f2py -c LB2D_np.f95 -m LB2Dw --f90flags="-O2 -static" --fcompiler=gnu95 --compiler=mingw32