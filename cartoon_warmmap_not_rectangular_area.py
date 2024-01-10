import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import math

#T=3010;
T=550
FFMpegWriter = manimation.writers['ffmpeg']
writer = FFMpegWriter(fps=25)
mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
writer.setup(fig, "video.mp4", 200)
#x = np.arange(0, np.pi*5, 0.1)
for i in range(1, T,10):
    print('Drawing frame %d...' % (i))
    with open("{id}.txt".format(id=i)) as file:
        data = [float(row.strip()) for row in file]

    z = data
    nrows = 3*40+1
    ncols = 2*40+1
    m=0
    subz=np.array(z)
    xlist = np.linspace(0.0,3.0,nrows)
    ylist = np.linspace(0.0,2.0,ncols)
    X, Y = np.meshgrid(xlist, ylist)
    Z = np.sqrt(X + Y)
    for k in range(0,ncols):
        for l in range(0,nrows):
            if (k>ncols/2 and l<nrows/3-1):
                Z[k,l]=math.nan
            else:
                Z[k,l]=subz[m]
                m+=1
    #Z = subz.reshape(ncols,nrows)
    #y = np.arange(ncols + 1)
    #x = np.arange(nrows + 1)
    plt.contourf(X,Y,Z,200)
    #plt.pcolormesh(x, y, Z, shading='flat', vmin=Z.min(), vmax=Z.max())
    #plt.subplots_adjust(bottom=0.1, right=0.95, top=0.9,wspace=1)
    plt.colorbar()
    #plt.plot(x, y)
    plt.title("T = %f" % (i))

    writer.grab_frame()

    plt.clf()