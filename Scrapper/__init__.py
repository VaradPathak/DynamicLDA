from lxml import html
import requests


page = requests.get('http://www.reuters.com/resources/archive/us/20080901.html')
tree = html.fromstring(page.text)

# This will create a list of article URLs:
URL = tree.xpath('//div[@class="headlineMed"]/a/@href')
Title = tree.xpath('//div[@class="headlineMed"]/a/text()')

print 'URLs: ', URL
print 'Title: ', Title
