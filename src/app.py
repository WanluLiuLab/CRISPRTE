from flask import Flask, render_template
import os

app = Flask(__name__)
app.secret_key = 'a04b1600-822a-11eb-97f9-acde48001122'