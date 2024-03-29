import gspread
from oauth2client.service_account import ServiceAccountCredentials
import os
import datetime
#Authorize the API

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
  sheet = client.open('Suggested Studies').sheet1
  entries = sheet.get_all_records()
  print('Sucessfully loaded suggested study file')
  if sheet.get('A1')[0][0] != 'Title':
    sheet.append_row(['Title', 'PMID', 'GEO ID', 'Contrast Conditions', 'Model System', 'Analysis Platform', 'Keywords', 'Authors', 'Contact Email'])
except:
  print('Error Authorizing Google Drive')



def log_suggested_study(user_data):
  try:
    row = user_data
    sheet.append_row(row)
    return "success"
  except Exception as e:
    print(e)
    print('Error adding to chat data to log')
    return 'failed'