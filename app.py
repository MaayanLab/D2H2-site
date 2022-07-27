from flask import Flask, render_template, url_for, request, jsonify
from pip import main
from query import *


app = Flask(__name__)


@app.route("/")
def index():
    return render_template("index.html")

@app.route('/getgwas', methods=['GET','POST'])
def get_gwas():

    gene = request.form['gene']

    if gene == '':
        return {'GWAS_Catalog':[]}

    result = query_gwas(gene)

    return result


@app.route('/getkomp', methods=['GET','POST'])
def get_mgi():

    gene = request.form['gene']
    result = query_mgi(gene)

    return result
