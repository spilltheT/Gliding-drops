import numpy as np
import os
import subprocess as sp
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import sys
import math

matplotlib.rcParams['text.usetex'] = True

def gettingFacets(filename, tracer):
	
    print('Getting facets values')
    if tracer == 1:
        exe = ["./getFacet1", filename]
    else:
        exe = ["./getFacet2", filename]

    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    segs = []
    skip = False
    if (len(temp2) > 1e1):
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                skip = False
                pass
            else:
                if not skip:
                    temp4 = temp2[n1+1].split(" ")
                    r1, z1 = np.array([float(temp3[1]), float(temp3[0])])
                    r2, z2 = np.array([float(temp4[1]), float(temp4[0])])                
                    segs.append(((r1, z1),(r2,z2)))
                    skip = True
    print('Got facets values for tracer')
    return segs

# ------------------------------------------------------------------------------

nGFS = 2500
lw = 3

angDiv = 6
slope = math.tan(math.pi/angDiv)

folder = 'Facets'  # output folder
if not os.path.isdir(folder):
    os.makedirs(folder)

rmin, rmax, zmin, zmax = [20, 38, -2, 3]

for ti in range(nGFS):
    t = 0.05 * ti
    place = "intermediate/snapshot-%5.4f" % t
    name = "%s/%8.8d.png" %(folder, int(t*1e3))
    if not os.path.exists(place):
        print("File %s not found!" % place)
    else:
        if os.path.exists(name):
            print("Image %s found!" % name)
        else:
            facets1 = gettingFacets(place, 1)
            if (len(facets1)):
                facets2 = gettingFacets(place, 2)

                fig, ax = plt.subplots()
                fig.set_size_inches(19.20, 10.80)
                rc('axes', linewidth=2)
                x_vals = np.linspace(0.0, rmax, num=20)
                y_vals = slope * x_vals
                ax.plot(x_vals, y_vals, '-',color='black',linewidth=lw/2)
                # ax.plot([0, rmin], [0, -rmin],'-',color='black',linewidth=lw/2)
                # ax.plot([0, rmax], [0, rmax],'-',color='black',linewidth=lw/2)
                ax.plot([rmin, rmin], [zmin, zmax],'-',color='black',linewidth=lw)
                ax.plot([rmin, rmax], [zmin, zmin],'-',color='black',linewidth=lw)
                ax.plot([rmin, rmax], [zmax, zmax],'-',color='black',linewidth=lw)
                ax.plot([rmax, rmax], [zmin, zmax],'-',color='black',linewidth=lw)

                line_segments2 = LineCollection(facets2, linewidths=4, colors='green', linestyle='solid')
                ax.add_collection(line_segments2)

                line_segments1 = LineCollection(facets1, linewidths=4, colors='orange', linestyle='solid')
                ax.add_collection(line_segments1)

                ax.set_aspect('equal')
                ax.set_xlim(rmin, rmax)
                ax.set_ylim(zmin, zmax)

                ax.set_title(r'$t/t_{\gamma} = %3.2f$' % t, fontsize=30)

                ax.axis('off')

                plt.savefig(name, bbox_inches="tight")
                plt.close()
            else:
                print("Problem in the available file %s" % place)

    print(("Done %d of %d" % (ti+1, nGFS)))
