import os
import pandas as pd
from dotenv import load_dotenv
import os
import tweepy
import re
import requests
 
# Use load_env to trace the path of .env:
load_dotenv('.env') 

auth = tweepy.OAuthHandler(os.environ.get("API_KEY"), os.environ.get("API_KEY_SECRET"))
auth.set_access_token(os.environ.get("ACCESS_TOKEN"), os.environ.get("ACCESS_TOKEN_SECRET"))

api = tweepy.API(auth)

try:
    api.verify_credentials()
    print("Authentication OK")
except:
    print("Error during authentication")

def update_tweets_table():
    tweets = api.user_timeline(screen_name="D2H2Bot", 
                            # 200 is the maximum allowed count
                            count=200,
                            include_rts = False,
                            tweet_mode = 'extended'
                            )

    hashtags_ignore = ['#D2H2', '#Enrichr', '#harmonizome', '#Harmonizome', "#T2D", '#diabetes']
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

    def parse_tweet(tweet_text) -> list:
        words = re.split("[\\n\\s]", tweet_text)
        seen_link = False
        seen_author = False
        prev_word = ''
        for i, word in enumerate(words):
            if '#' in word and not(word in hashtags_ignore) and 'http' not in word:
                gene = word.replace('#', '')
            elif '@' in word:
                journal = journal_dict[word]
            elif prev_word == 'by' and not seen_author:
                author = " ".join(words[i:i+3])
                title_start = i+3
                seen_author = True
            elif 'https' in word and not(seen_link):
                title = " ".join(words[title_start:i]).replace('\n', '')
                article_link = words[i]
                seen_link = True  
            prev_word = word
        article_link = requests.head(article_link).headers['location']
        enrichr = 'https://maayanlab.cloud/Enrichr/#find!gene=' + gene
        harmonizome = 'https://maayanlab.cloud/Harmonizome/gene/' + gene
        D2H2 = 'https://d2h2.dev.maayanlab.cloud/singlegene/' + gene
        return [gene, title, author, journal, article_link, [D2H2, enrichr, harmonizome]]


    tweets_table = []
    for info in tweets:
        text = info.full_text
        tweet_list = parse_tweet(text)
        tweets_table.append(tweet_list)


    df = pd.DataFrame(tweets_table, columns=["Gene","Title","Author(s)", "Journal" ,"Article", "Analyze"])
    df.to_csv('static/searchdata/tweets.csv', index=False)











