import os
import sys


def buildtriple(idprob,vocab):
    (word_id, prob) = idprob.split(':')
    word = vocab[int(word_id)]
    return (word, word_id, prob)
    
#input: directory of input files, directory for output files, vocab file
if len(sys.argv) >= 3:
    infolder = sys.argv[1]
    outfolder = sys.argv[2]
    vocabfilename = sys.argv[3]

    vocabfile = open(vocabfilename, 'r')
    vocab = [word.strip() for word in vocabfile.readlines()]
    
    #iterate through files of infolder
    for filename in os.listdir(infolder):
        thefile = open(infolder + '/' + filename)
        for topic in thefile.readlines():
            outline = [buildtriple(idprob,vocab) for idprob in topic.split(',')]
            map(lambda x: thefile.write(str(x)), outline)
            thefile.write('\n')
