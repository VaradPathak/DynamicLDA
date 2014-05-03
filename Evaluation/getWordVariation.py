import os
import sys

if len(sys.argv) >= 3:
    topic_id = int(sys.argv[1])
    word_id = int(sys.argv[2])
    filepath = sys.argv[3]
    startYear = int(sys.argv[4])
    endYear = int(sys.argv[5])

files = os.listdir(filepath)

for year in range(startYear, endYear + 1):
    for month in range(1, 13):
        if(month > 9):
            fileName = filepath + '/' + 'topics_' + str(year) + str(month) + '.txt'
            if os.path.isfile(fileName):
#                 print 'reading file: ' + fileName
                monthFile = open(fileName, 'r')
                lines = monthFile.readlines()
                topic = lines[topic_id]
                prob = float(topic.split(',')[word_id])
                print prob
        else:
            fileName = filepath + '/' + 'topics_' + str(year) + '0' + str(month) + '.txt'
            if os.path.isfile(fileName):
#                 print 'reading file: ' + fileName
                monthFile = open(fileName, 'r')
                lines = monthFile.readlines()
                topic = lines[topic_id]
                prob = float(topic.split(',')[word_id])
                print prob
