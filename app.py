from flask import Flask, request, Response
import json
from flask_cors import CORS, cross_origin
from maria_model import predict_with4, make_mhc_list, scan_gene

app = Flask(__name__)
CORS(app)

with open("index.html") as f:
    html_page = f.read()

with open("bundle.js") as f:
    js_bundle = f.read()

@app.route("/")
def index():
    return html_page

@app.route("/bundle.js")
def js():
    with open("bundle.js") as f:
        js_bundle = f.read()
    return js_bundle

@app.route("/predict", methods=["POST"])
def predict():
    #potential error message
    error0 = []

    #get data
    data = json.loads(request.data)
    data_dr = data["data"][0]
    data_gn = data["data"][1]
    data_seq = data["data"][2]

    #process data
    pep_list = [pep0.rstrip() for pep0 in data_seq.split("\n")]
    if len(pep_list[-1]) == 0:
        pep_list = pep_list[0:-1]
    gn_list = [gn0.rstrip() for gn0 in data_gn.split("\n")]
    if len(gn_list[-1]) == 0:
        gn_list = gn_list[0:-1]

    #scan gene when the input is a single gene without the peptide info
    if len(pep_list) == 0 and len(gn_list) == 1:
        list_out = scan_gene(gn_list[0])
        if len(list_out) > 1:
            gn_list = list_out[0]
            pep_list = list_out[1]
        else:
            error0.append(list_out)

    #check error 
    if not len(pep_list) == len(gn_list):
        error0.append('Error: Inconsistent gene and peptide numbers.')
    mhc_list = make_mhc_list(data_dr.rstrip(),pep_list)
    if 'allele' in mhc_list:
        error0.append(mhc_list)

    #only run MARIA if the input data is valid 
    if len(error0) > 0:
        output = error0
    else:
        #run maria 
        output = pep_list + [str(round(x,4)) for x in list(predict_with4(pep_list,mhc_list,gn_list,ranking=False))]
        #print(output)
        #print(error0)

    return json.dumps(output)
    #[len(x.split()) for x in data["data"].split("\n")]

if __name__ == "__main__":
    app.run()

