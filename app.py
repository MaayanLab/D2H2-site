from flask import Flask, render_template, url_for, request, jsonify
from pip import main
from query import *


app = Flask(__name__)


@app.route("/")
def index():
    return render_template("index.html")

@app.route('/getgwas', methods=['GET','POST'])
def get_enricher():

    gene = request.form['gene']

    if gene == '':
        return {'GWAS_Catalog_2019':[]}

    result = query_enricher(gene)

    return result


@app.route('/getkomp', methods=['GET','POST'])
def get_komp():

    gene = request.form['gene']
    result = query_komp(gene)

    return result
