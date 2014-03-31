from lxml import html
import requests
import article as ac


docId = 0  
f = open('2014.txt', 'w')
for yr in range(2013, 2014):
    year = 'http://www.reuters.com/resources/archive/us/' + str(yr)
    for mnth in range(1, 2):
        if(mnth < 10):
            month = '0' + str(mnth)
        else:
            month = str(mnth)
        
        for day in range(1, 2):
            if(day < 10):
                URL = year + month + '0' + str(day) + '.html'
            else:
                URL = year + month + str(day) + '.html'
            
            page = requests.get(URL)
            tree = html.fromstring(page.text)
        # This will create a list of article URLs:
            URLs = tree.xpath('//div[@class="headlineMed"]/a/@href')
        #     Title = tree.xpath('//div[@class="headlineMed"]/a/text()')
            date = URL[-13:-5]
            for num in range(0, len(URLs)):
                docId = docId + 1
                doc = ac.article('', date, '', URLs[num], docId)
                curpage = requests.get(doc.URL)
                curtree = html.fromstring(curpage.text)
                Title = curtree.xpath('//*[@id="content"]/div[4]/div/div[3]/div[1]/h1/text()')
                Paragraphs = curtree.xpath('//*[@id="articleText"]/p/text()')
                # Location = tree.xpath('//*[@id="articleInfo"]/p[2]/span[1]/text()')
                if len(Title) > 0:
                    doc.Title = Title[0].replace('\"', '')
                    Paragraphs.append(Title[0])
                doc.Text = " ".join(Paragraphs)
                doc.Text = doc.Text.replace('\n', ' ')
                doc.Text = doc.Text.replace('\"', '')
                f.write(str(doc.id) + "," + str(doc.Date) + ',\"' + doc.Title.encode('utf-8') + '\",\"' + doc.Text.encode('utf-8') + '\"\n')

f.close()
    
