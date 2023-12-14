import math
import numpy as np
from params import nz, step, istart, iend, ignore_before

yy=np.zeros(nz)
avg=np.zeros(nz)
avg_tot=np.zeros(nz)
avg_stress=np.zeros(nz)
avg_stress1=np.zeros(nz)
avg_stress2=np.zeros(nz)
avg_urms=np.zeros(nz)
avg_vrms=np.zeros(nz)
avg_wrms=np.zeros(nz)

count=0

with open("profiles-2.dat") as f:
    for istep in range(istart, iend+step, step):
        line=f.readline()
        fstep=int(line.strip().split()[0])

        if istep>ignore_before:
                count=count+1
        
        for k in range(0, nz):
            line=f.readline()
            y0=k+1
            row=line.strip().split()
            ux=float(row[0])
            uy=float(row[1])
            uz=float(row[2])

            stress_1=float(row[3])
            stress_2=float(row[4])
            stress_3=float(row[5])

            urms=float(row[6])
            vrms=float(row[7])
            wrms=float(row[8])

            yy[k]=y0
            
            if istep>ignore_before:
                avg[k]=avg[k]+uz
                avg_stress[k]=avg_stress[k]+stress_1
                avg_urms[k]=avg_urms[k]+urms
                avg_vrms[k]=avg_vrms[k]+vrms
                avg_wrms[k]=avg_wrms[k]+wrms


    avg=avg/(count)
    avg_stress=avg_stress/count
    avg_urms=avg_urms/count
    avg_vrms=avg_vrms/count
    avg_wrms=avg_wrms/count

    # now merge two halfs
    avg_rev=avg[::-1]
    avg_tot=avg
    avg_tot_rev=avg[::-1]
    avg=0.5*(avg+avg_rev)

    avg_rev=avg_stress[::-1]
    avg_stress1=avg_stress
    avg_stress2=-avg_stress[::-1]
    avg_stress=0.5*(avg_stress-avg_rev)

    avg_rev=avg_urms[::-1]
    avg_urms=0.5*(avg_urms+avg_rev)
    avg_rev=avg_vrms[::-1]
    avg_vrms=0.5*(avg_vrms+avg_rev)
    avg_rev=avg_wrms[::-1]
    avg_wrms=0.5*(avg_wrms+avg_rev)

    # save mean velocity
    with open("VELMEAN.dat", "w") as g:
        for k in range(int(nz)+0):
            print(k)
            g.write(f"{yy[k]} {avg[k]} {avg_stress[k]} {avg_urms[k]} {avg_vrms[k]} {avg_wrms[k]}\n")
