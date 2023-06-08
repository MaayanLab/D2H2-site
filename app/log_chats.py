import gspread
from oauth2client.service_account import ServiceAccountCredentials
import dotenv
import os
import datetime

dotenv.load_dotenv()
#Authorize the API

DEBUG = os.environ.get('DEBUG', True)
scope = [
    'https://www.googleapis.com/auth/drive',
    'https://www.googleapis.com/auth/drive.file'
    ]
file_name = 'client_key.json'

creds_dict = {
  "type": "service_account",
  "project_id": os.environ.get('PROJECT_ID'),
  "client_email": os.environ.get('CLIENT_EMAIL'),
  "client_id": os.environ.get('CLIENT_ID'),
  "auth_uri": "https://accounts.google.com/o/oauth2/auth",
  "token_uri": "https://oauth2.googleapis.com/token",
  "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
  "client_x509_cert_url": os.environ.get('CERT_URL'),
  "universe_domain": "googleapis.com",
  "private_key_id": os.environ.get('PRIVATE_KEY_ID'),
  "private_key": os.environ.get('PRIVATE_KEY').replace('\\n', '\n')
}

try:
  creds = ServiceAccountCredentials.from_json_keyfile_dict(creds_dict,scope)
  client = gspread.authorize(creds)
  sheet = client.open('Chat logs').sheet1
  entries = sheet.get_all_records()
  print('Sucessfully loaded chat logs')
except:
  print('Error Authorizing Google Drive')



def log_chat(user_query, response):
  if not DEBUG:
    try:
      row = [str(datetime.datetime.now()), user_query, response]
      sheet.append_row(row)
      return
    except:
      print('Error adding to chat data to log')
      return
