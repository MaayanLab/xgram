
from flask import Flask
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash
import json
import sys
import logging
from logging.handlers import RotatingFileHandler
import os
from flask import send_from_directory

# # change routing of logs when running docker 
# logging.basicConfig(stream=sys.stderr) 



# app = Flask(__name__)
app = Flask(__name__, static_url_path='')

ENTRY_POINT = '/chikungunyaKEA'

# switch for local and docker development 
# docker_vs_local
##########################################

# for local development 
SERVER_ROOT = os.path.dirname(os.getcwd()) + '/chikungunyaKEA/chikungunyaKEA' ## original 

# # for docker development
# SERVER_ROOT = '/app/chikungunyaKEA'


@app.route(ENTRY_POINT + '/<path:path>') ## original 
# @crossdomain(origin='*')
def send_static(path):

  # print('path and SERVER_ROOT')
  # print(path)
  # print(SERVER_ROOT)
  
  return send_from_directory(SERVER_ROOT, path)


@app.route("/chikungunyaKEA/")
def index():
  print('Rendering index template')
  return render_template("index.html")

# post request
@app.route('/chikungunyaKEA/', methods=['GET','POST'])
def python_function():
  import flask 
  # import make_enr_clust
  import make_exp_clustergram
  import json 
  import json_scripts 
  import sys
  import cookielib, poster, urllib2, json

  error = None 

  # get request json 
  req_json = request.form

  # get the gene class 
  inst_name = req_json['name'] # request.form['name'].upper()

  # generate an expression clustergram using cst code 
  d3_json = make_exp_clustergram.tf_clust(inst_name)

  # jsonify network 
  #####################
  return flask.jsonify(d3_json)

if __name__ == "__main__":
    app.run(host='0.0.0.0',port=5000,debug=True)
 