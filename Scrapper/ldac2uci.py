infile = open('jan.dat')
outfile = open('JanUCI.txt', 'w')

doclines = infile.readlines()
infile.close()

for i in range(len(doclines)):
    line = doclines[i].strip().split(' ')[1:]
    for elt in line:
        pieces = elt.split(':')
        outfile.write(str(i + 1) + ' ' + str(int(pieces[0]) + 1) + ' ' + pieces[1] + '\n')
    
outfile.close()
