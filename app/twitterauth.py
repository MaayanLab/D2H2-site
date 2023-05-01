import pandas as pd
import os
import re
import requests
import xml.etree.ElementTree as ET
from functools import lru_cache




def parse_tweet(tweet_text, html) -> list:
    hashtags_ignore = ['#D2H2', '#Enrichr', '#harmonizome', '#Harmonizome', "#T2D", '#diabetes', '(@D2H2Bot)']
    journal_dict = {
    '@Nature': 'Nature', 
    '@Cell_Metabolism': 'Cell Metabolism', 
    '@ScienceMagazine': 'Science',
    '@PNASNews': 'PNAS',
    '@NEJM': 'NEJM',
    '@JAMA_current': 'JAMA',
    '@NatureMedicine': 'Nature Medicine',
    '@NatureGenet': 'Nature Genetics',
    '@NatureCellBio': 'Nature Cell Biology',
    '@NatureBiotech': 'Nature Biotechnology',
    '@NatMetabolism': 'Nature Metabolism',
    '@CellCellPress': 'Cell'
    }
    words = re.split("[\\n\\s]", tweet_text)
    if words[-1].count('#') > 1:
        hashtags = words[-1]
        hashtags = hashtags.split('#')[1:]
        words.remove(words[-1])
        words = words + list(map(lambda x: '#' + x, hashtags))
    title_done = False
    seen_author = False
    prev_word = ''
    for i, word in enumerate(words):

        if '#' in word and not(word in hashtags_ignore) and 'http' not in word:
            gene = word.replace('#', '')
            title_range = words[title_start:i]
            title = " ".join(filter(lambda x: 'http' not in x and '#' not in x, title_range)).replace('\n', '')
            article_link= list(filter(lambda x: 'http' in x, title_range))[0]
        elif '@' in word and word not in hashtags_ignore:
            journal = journal_dict[word]
        elif prev_word == 'by' and not seen_author:
            author = " ".join(words[i:i+3])
            title_start = i+3
            seen_author = True
            
        prev_word = word
    title = title.strip()

    #article_link = html.split(title.split(' ')[-1] + '<')[1].split('href=')[1].split('\"')[1]
    #article_link = requests.head(article_link).headers['location']
    enrichr = 'https://maayanlab.cloud/Enrichr/#find!gene=' + gene
    harmonizome = 'https://maayanlab.cloud/Harmonizome/gene/' + gene
    D2H2 = 'https://d2h2.dev.maayanlab.cloud/singlegene/' + gene
    return [gene, title, author, journal, article_link, [D2H2, enrichr, harmonizome]]


@lru_cache()
def update_tweets_table(day):
    tweets_table = []
    try:
        xml_text = requests.get('https://nitter.net/D2H2Bot/rss', headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_5) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/50.0.2661.102 Safari/537.36'}).text
        root = ET.fromstring(xml_text)
    except:
        print('error fetching twitter data')
        return

    tweets = []
    for item in root[0].findall('item'):
        t =[]
        t.append(item.find('title').text)
        t.append(item.find('description').text)
        tweets.append(t)

    for tweet in tweets:
        text = tweet[0]
        html = tweet[1]
        tweet_list = parse_tweet(text, html)
        tweets_table.append(tweet_list)

    df = pd.DataFrame(tweets_table, columns=["Gene","Title","Author(s)", "Journal" ,"Article", "Analyze"])
    df.to_csv('static/searchdata/tweets.csv', index=False)
