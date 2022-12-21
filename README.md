D2H2 is a Flask application that has both Python and R dependcies. Please first ensure you have Python 3.9 installed and if not please download through homebrew or the offical python site.

Create a Python 3.9 enviroment and install the dependencies outlined in requirements.txt.

```
python3.9 -m venv env
source env/bin/activate
pip3 install -r requirements.txt
```

For the differential gene analysis, R is used through the module rpy2. Please ensure R is installed and run setup.R script contained within the app directory to ensure the relvevant depenencies are installed (Bioconductor, limma, edgeR, DESeq2, R.utils, RCurl, statmod).

You must also run the redis message broker server. This will enable tasks to be sent to the celery worker running in the background. If you don't already have redis installed you can install and run it with the following commands:

```
wget http://download.redis.io/redis-stable.tar.gz
tar xvzf redis-stable.tar.gz
cd redis-stable
make
sudo make install
```

Then you can run the redis server from the redis-stable directory:

```
redis-server
```

To start the celery worker that will excecute the differential gene expression analysis from the D2H2-site directory (make sure the enviroment with the dependencies is activated):

```
cd app
celery -A celery_config.celery worker --loglevel=info -P eventlet
```

Then to run the application, in a separate terminal, navigate to the the app directory and run the file app.py:

```
cd app
python3 app.py
```


