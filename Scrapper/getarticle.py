# -*- coding: utf-8 -*-
from lxml import html
import requests
#http://www.reuters.com/article/2012/12/02/us-space-france-russia-idUSBRE8B101L20121202
#http://www.reuters.com/article/2007/01/02/music-jazz-chicago-dc-idUSN2927338620070102
#http://www.reuters.com/article/2014/03/28/us-microsoft-office-ipad-idUSBREA2Q1MV20140328
page = requests.get('http://www.reuters.com/article/2014/01/02/walmart-china-idUSL3N0KC0LH20140102')
tree = html.fromstring(page.text)

# This will create a list of article URLs:
#URL = tree.xpath('//div[@class="headlineMed"]/a/@href')
#Title = tree.xpath('//div[@class="headlineMed"]/a/text()'
Title = tree.xpath('//*[@id="content"]/div[4]/div/div[3]/div[1]/h1/text()')

Location = tree.xpath('//*[@id="articleInfo"]/p[2]/span[1]/text()')

Paragraphs = tree.xpath('//*[@id="articleText"]/p/text()')

print 'Paragraphs: ', Paragraphs
print 'Location: ' , Location
print 'Title:' , Title
