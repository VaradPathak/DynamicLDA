import os


vocabDict = {}

files = os.listdir('/home/vspathak/git/DynamicLDA/Scrapper/output')

for fileName in files:
    for line in open('output/'+fileName, 'r'):
        line = line.strip()
        words = line.split()
        for key in words:
            if key in vocabDict:
                vocabDict[key] += 1
            else:
                vocabDict[key] = 1

stopwordsFile = open('stopwords', 'a')
for term in vocabDict:
    if vocabDict[term] < 26:
        print term
        stopwordsFile.write(term.encode('utf-8') + '\n')

stopwordsFile.close()
