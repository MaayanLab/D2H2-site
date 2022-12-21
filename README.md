# D2H2

D2H2 is a Flask application that has both Python and R dependcies. Please first ensure you have Python 3.9 installed and if not please download through homebrew or the offical python site. 

Create a Python 3.9 enviroment and install the dependencies outlined in requirements.txt.

''' 
python3.9 -m venv env
source env/bin/activate
pip3 install -r requirements.txt
'''

For the differential gene analysis, R is used through the module rpy2. Please ensure R is installed and run setup.R script contained within the app directory to ensure the relvevant depenencies are installed (Bioconductor, limma, edgeR, DESeq2, R.utils, RCurl, statmod).

Then to run the application, navigate to the the app directory and run the file app.py.

'''
cd app
python3 app.py
'''





