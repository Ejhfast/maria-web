USEFUL COMMANDS:
npm install # install everything in package.json
npm install PACKAGE --save # adds to the package.json
webpack # compile the javascript in entry.js into bundle.js

FILES:
bundle.js: compiled bundle, don't ever edit
entry.js: javascript views/app
index.html: dummy page, don't need to edit
app.py: flask app that serves JS and provides Maria APIs
package.json: contains npm/JS dependencies for project
.babelrc: babel setting, allows use of ES6 and JSX

SUMMARY:
Only ever need to edit app.py and entry.js
