import os
import sys


def buildtriple(idprob,vocab):
#     if (len(idprob) == 0):
#         print 'working'
    word_Prob = idprob.split(':')
    word = vocab[int(word_Prob[0])]
    return (word, word_Prob[0], word_Prob[1])
    
#input: directory of input files, directory for output files, vocab file
if len(sys.argv) >= 3:
    infolder = sys.argv[1]
    outfolder = sys.argv[2]
    vocabfilename = sys.argv[3]

    vocabfile = open(vocabfilename, 'r')
    vocab = [word.strip() for word in vocabfile.readlines()]
    vocabfile.close()
    
    #iterate through files of infolder
    for filename in os.listdir(infolder):
        infile = open(infolder + '/' + filename)
        outfile = open(outfolder + '/' + filename, 'w')
        for topic in infile.readlines():
            temp = topic[:-2].split(',')
            word_Prob = temp[0].split(':')
            if not (word_Prob[0].isdigit()):
                continue
            outline = [buildtriple(idprob,vocab) for idprob in topic[:-2].split(',')]
            map(lambda x: outfile.write(str(x)), outline)
            outfile.write('\n')
        outfile.close()
        infile.close()
