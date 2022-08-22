import requests
from bs4 import BeautifulSoup
import re

URL = "https://pubmed.ncbi.nlm.nih.gov/19304878"
page = requests.get(URL)

soup = BeautifulSoup(page.content, "html.parser")


print(soup)

results = soup.findAll("title")

print(results[0].text)

