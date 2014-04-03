infile = open('reuters1monthsub/reuters1monthsub.dat')
outfile = open('reuters1monthsubUCI.csv','w')

doclines = infile.readlines()
infile.close()

for i in range(len(doclines)):
    line = doclines[i].strip().split(' ')[1:]
    for elt in line:
        pieces = elt.split(':')
        outfile.write(str(i+1) + ' ' + pieces[0] + ' ' + pieces[1] + '\n')
    
outfile.close()



