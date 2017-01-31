from flask import Flask, request, Response
import json

app = Flask(__name__)

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
    data = json.loads(request.data)
    print(data)
    return json.dumps([len(x.split()) for x in data["data"].split("\n")])

if __name__ == "__main__":
    app.run()
