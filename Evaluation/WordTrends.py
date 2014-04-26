import os

topTenFile = open('/home/vspathak/git/DynamicLDA/TopTen.txt', 'r')
for input in topTenFile:
    words = input.split(',')
    topSet = set(words)

wordEvolution = open('/home/vspathak/git/DynamicLDA/WordEvolution.txt', 'w')

files = os.listdir('/home/vspathak/git/DynamicLDA/Scrapper/Pi')

for fileName in files:
    for topic in open('output/' + fileName, 'r'):
        topic = topic.strip()
        wordProb = topic.split(',')
        wordEvolution.write(wordProb[word])
