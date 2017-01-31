import React from 'react';
import { render } from 'react-dom';
import fetch from 'isomorphic-fetch';
import axios from 'axios';

const postData = (input, component) => {
  axios.post('http://localhost:5000/predict', {
    data: input.value
  })
  .then(response => {
    return response.data;
  })
  .then(json => {
    component.setState({
      scores: json
    })
  })
};

let input;

class DataEntry extends React.Component {
  constructor (props) {
    super(props)
    this.state = {scores: []};
  }
  render () {
    const styles = this.constructor.styles;
    return (<div>
        <div className="left_pane" style={styles.pane}>
          <textarea className="input_data" style={styles.textarea} ref={node => {input = node;}}></textarea>
          <br />
          <button style={styles.button} onClick={() => postData(input, this)}>Analyze</button>
        </div>
        <DataDisplay scores={this.state.scores}/>
      </div>);
  }
}

DataEntry.styles = {
  pane: {
    margin: "20px 0px"
  },
  textarea: {
    width: "95%",
    height: "200px",
    margin: "0px 0px 10px 0px",
    padding: "1%",
    border: "1px solid #ccc",
    fontSize: ".9em"
  },
  button: {
    fontSize: "1em"
  }
}

class DataDisplay extends React.Component {
  render () {
    const styles = this.constructor.styles;
    const pane_style = this.props.scores.length == 0 ? { display: "none" } : styles.pane;
    const format_scores = this.props.scores.join("\n");
    return (<textarea className="right_pane" style={pane_style} value={format_scores}></textarea>);
  }
}

DataDisplay.styles = {
  pane: {
    width: "95%",
    height: "200px",
    padding: "1%",
    lineHeight: "1.5em",
    fontSize: ".9em",
    backgroundColor: "#444",
    fontFamily: "monaco, courier new",
    color: "#eee"
  }
}

class App extends React.Component {
  render () {
    const styles = this.constructor.styles;
    return (<div style={styles.body}>
        <h1 style={styles.h1}>Maria</h1>
        <h2 style={styles.h2}>A deep learning based tool for predicting peptide presentation</h2>
        <DataEntry />
      </div>);
  }
}

App.styles = {
  body: {
    fontSize: "14px",
    fontFamily: "helvetica neue, helvetica, sans-serif",
    width: "800px",
    margin: "20px auto"
  },
  h1: {
    fontSize: "1.5em",
    fontWeight: "800",
    margin: "3px 0px",
  },
  h2: {
    fontSize: "1.3em",
    fontWeight: "200",
    margin: "3px 0px"
  }
}

render(<App/>, document.getElementById('container'));
