from lxml import html
import requests
import article as ac


urls = list()
articles = list()
for yr in range(2007, 2008):
    year = 'http://www.reuters.com/resources/archive/us/' + str(yr)
    for mnth in range(1, 13):
        if(mnth < 10):
            month = '0' + str(mnth)
        else:
            month = str(mnth)
        for day in range(1, 32):
            if(day < 10):
                URL = year + month + '0' + str(day) + '.html'
            else:
                URL = year + month + str(day) + '.html'
            urls.append(URL)
    
docId = 0        
for URL in urls:
    page = requests.get(URL)
    tree = html.fromstring(page.text)
    
    # This will create a list of article URLs:
    URLs = tree.xpath('//div[@class="headlineMed"]/a/@href')
    Title = tree.xpath('//div[@class="headlineMed"]/a/text()')
    date = URL[-13:-5]
    
    for num in range (0, len(URLs)):
        docId = docId + 1 
        doc = ac.article(Title[num], date, '', URLs[num], docId)
        articles.append(doc)

for doc in articles:
    print doc.id, doc.Date, doc.URL, doc.Title
