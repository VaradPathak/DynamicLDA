import os
import sys

if len(sys.argv) >= 2:
    perpfilename = sys.argv[1]
    num_iter = int(sys.argv[2])

perpfile = open(perpfilename)

finalperps = [perp for perp in perpfile.readlines() if perp[1] != ',' and int(perp[:2]) == num_iter]

# for perp in finalperps:
#     print perp.split(',')[1]

#print len(finalperps)
for i in range(len(finalperps)):
    if i % 3 == 0:
        print i/3, finalperps[i].split(',')[1]
