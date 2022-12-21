import os
import s3fs
from celery import Celery 
from flask import Flask

endpoint = os.environ.get('ENDPOINT', 'https://minio.dev.maayanlab.cloud/')
base_url = os.environ.get('BASE_URL', 'd2h2/data')
ROOT_PATH = os.environ.get('ROOT_PATH', '/')
BASE_PATH = os.environ.get('BASE_PATH', 'maayanlab.cloud')
DEBUG = os.environ.get('DEBUG', True)

print(endpoint)
s3 = s3fs.S3FileSystem(anon=True, client_kwargs={'endpoint_url': endpoint})


app = Flask(__name__, static_url_path=ROOT_PATH + 'static')
print(app.name)

app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'

app.app_context().push()

celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'], include=["dge"])
celery.conf.update(app.config)