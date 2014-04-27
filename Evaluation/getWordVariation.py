import os
import sys

if len(sys.argv) >= 3:
    topic_id = int(sys.argv[1])
    word_id = int(sys.argv[2])
    filepath = sys.argv[3]

files = os.listdir(filepath)

for fileName in files:
    monthFile = open(filepath + '/' + fileName, 'r')
    lines = monthFile.readlines()
    topic = lines[topic_id]
    prob = float(topic.split(',')[word_id])
    print prob
