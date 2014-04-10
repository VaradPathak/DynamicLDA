from lxml import html
import requests
import article as ac
import sys
import random
import nltk
import re
from porter2 import stem


docId = 0  
if len(sys.argv) >= 4:
    theyear = int(sys.argv[1])
    firstmonth = int(sys.argv[2])
    num_months = int(sys.argv[3])
    seqfile = open('seq-' + str(theyear) + '-' + str(firstmonth) + '-' + str(num_months) + '.txt', 'w')

else:
    print 'usage: python __init__.py year firstmonth num_months'

for yr in range(theyear, theyear + 1):
    year = 'http://www.reuters.com/resources/archive/us/' + str(yr)
    for mnth in range(firstmonth, firstmonth + num_months):
        if(mnth < 10):
            month = '0' + str(mnth)
        else:
            month = str(mnth)
        
        monthdocs = 0
        for day in range(1, 32):
            if(day < 10):
                URL = year + month + '0' + str(day) + '.html'
            else:
                URL = year + month + str(day) + '.html'
            
            page = requests.get(URL)
            tree = html.fromstring(page.text)
            URLs = tree.xpath('//div[@class="headlineMed"]/a/@href')
            date = URL[-13:-5]

            f = open('output/' + str(date) + '.txt', 'w')
            # generate the random vector(python generate a sample without 
            # replacement from a range of numbers)
             
            for num in random.sample(range(0, len(URLs)), int(len(URLs) / 10)):
                doc = ac.article('', date, '', URLs[num], -1)
                curpage = requests.get(doc.URL)
                curtree = html.fromstring(curpage.text)
                Title = curtree.xpath('//*[@id="content"]/div[4]/div/div[3]/div[1]/h1/text()')
                Paragraphs = curtree.xpath('//*[@id="articleText"]/p/text()')
                if len(Title) > 0:
                    doc.Title = Title[0].replace('\"', '')
                    Paragraphs.append(Title[0])
                doc.Text = " ".join(Paragraphs)
                doc.Text = doc.Text.replace('\n', ' ')
                doc.Text = doc.Text.replace('\"', '')

                if(len(doc.Text.split()) > 100):
                    docId = docId + 1
                    doc.id = docId
                    print doc.id
                    monthdocs = monthdocs + 1
                    
                    docText = re.sub('[^A-Za-z]+', ' ', doc.Text)
                    docTitle = re.sub('[^A-Za-z]+', ' ', doc.Title)
                    docText = docTitle + ' ' + docText  
                    docText = docText.lower()
                    tokens = docText.split()

                    docText = " ".join([stem(t) for t in tokens])

                    f.write(docText.encode('utf-8') + '\n')

            f.close()
        seqfile.write(str(theyear) + '-' + str(mnth) + ':' + str(monthdocs) + '\n')
seqfile.close()
