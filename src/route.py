from flask import render_template, request, flash, redirect, url_for, session
from app import app
import random
import os
from os.path import join
import pandas as pd
import json
from utils import *
from grna_mismatch_scoring import *
import uuid
import multiprocessing

tb_name = "./database/Homo_sapiens.GRCh38.97.dna.primary_assembly.gRNA.db"
app.config['RESULTS_FOLDER'] = './results'
app.config['cpu_count'] = multiprocessing.cpu_count()

@app.route('/crisprte', methods=["GET"])
def index():
    return render_template("index.html")


@app.route('/crisprte/documention', methods=["GET"])
def instructions():
    return render_template("documentation.html")


@app.route('/crisprte/result', methods=["GET"])
def result():
    form_data = request.values.to_dict()
    if 'uid' in form_data.keys() and form_data['uid'] != '':
        session['uid'] = form_data["uid"]
        if session['uid'].startswith("A"):
            return render_template("result.html", uid=form_data["uid"])
        elif session['uid'].startswith("B"):
            return render_template("result_combination.html", uid=form_data["uid"])
    elif form_data['for'] == 'K2':
        uid = uuid.uuid4().hex + "A"
        session['uid'] = uid
        return render_template("result.html", TE_dup=form_data["te"], message="We only shows 100 gRNAs for knockouting TE subclass. Please query knockout single duplicate for more gRNAs.", uid=uid)
    elif form_data['for'] == 'K3':
        uid = uuid.uuid4().hex + "B"
        session['uid'] = uid
        return render_template("result_combination.html", TE_dup=form_data["te"], message="The combination aims to cover maximum copies from a TE subclass based on random sampling.", uid=uid)
    elif form_data['for'] == 'K1':
        uid = uuid.uuid4().hex + "A"
        session['uid'] = uid
        return render_template("result.html", TE_dup=form_data["te"], uid=uid)


@app.route('/crisprte/api', methods=["POST"])
def api():
    if request.method == "POST":
        form_data = request.json
        # print(form_data)
        result = None
        if form_data['type'] == 'gtf':
            conn = SQLite.connect("./database/hg38.db")
            cursor = conn.cursor()
            result = SQLite.select_gtf_by_region(cursor, 'pcg', form_data['chrom'], form_data['start'], form_data['end']) + SQLite.select_te_gtf_by_region(
                cursor, 'te', form_data['chrom'], form_data['start'], form_data['end'])
            conn.close()
            result = json.dumps(result)

        if form_data['type'] == 'uid':
            with open(join(app.config['RESULTS_FOLDER'], form_data['uid'] + ".json"), "r") as f:
                result = json.load(f)
                return result

        elif form_data['type'] == 'te':
            with open("./database/hg38_TEs.json") as f:
                return json.load(f)
            

        elif form_data['type'] == 'K1':
            mms = query_mismatch_annotation(tb_name, te_dup=form_data['te'])
            result = calculate_offtarget_score(
                mms, form_data['te'].split("_dup")[0], form_data['te'])
            with open(join(app.config['RESULTS_FOLDER'], session['uid'] + ".json"), "w+") as f:
                f.write(json.dumps({"data":result,"te": form_data['te']}))
            result = json.dumps(result)


        elif form_data['type'] == 'K2':
            mms = query_mismatch_annotation(tb_name, te_class=form_data['te'])
            result = calculate_offtarget_score(
                mms, form_data['te'].split("_dup")[0], form_data['te'])
            with open(join(app.config['RESULTS_FOLDER'], session['uid'] + ".json"), "w+") as f:
                f.write(json.dumps({"data":result,"te": form_data['te']}))
            result = json.dumps(result)


        elif form_data['type'] == 'K3':
            gids_set = None
            coverage = None
            result = {}
            combination = get_combination(tb_name, gids_set)
            calculate_offtarget_score(combination, form_data['te'], None)
            result["data"] =  {"data": combination, "coverage": coverage, "gids_set": gids_set}
            result["te"] = form_data['te']

            result = json.dumps(result)
            with open(join(app.config['RESULTS_FOLDER'], session['uid'] + ".json"), "w+") as f:
                f.write(result)
        return result
