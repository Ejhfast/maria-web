import React from 'react';
import { render } from 'react-dom';
import fetch from 'isomorphic-fetch';
import axios from 'axios';

const postData = (input_allele, input_gene, input_seq, component ) => {
  axios.post('/predict', {
    data: [input_allele.value, input_gene.value, input_seq.value, ]
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

let input_allele,input_gene,input_seq; 

class DataEntry extends React.Component {
  constructor (props) {
    super(props)
    this.state = {scores: []};
  }
  render () {
    const styles = this.constructor.styles;
    return (<div>
        <div className="allele_pane" style={styles.pane}>   
          <textarea className="input_data" style={styles.textarea} ref={node => {input_allele = node;}}></textarea>
          <div className="label" style={styles.subtitle} > <b>HLA-DR allele: </b> Follow the format HLA-DRB1*07:01, single or double alleles (separtaed by a comma).</div>
        </div>
        <div className="gene_pane" style={styles.pane}>
          <textarea className="input_data" style={styles.textarea} ref={node => {input_gene = node;}}></textarea>
          <div className="label" style={styles.subtitle} >  <b>Gene: </b> HUGO gene symblos for the corresponding gene of peptides (e.g. KRAS, one per line).</div>
        </div>
        <div className="seq_pane" style={styles.pane}>
          <textarea className="input_data" style={styles.textarea} ref={node => {input_seq = node;}}></textarea>
        <div className="label" style={styles.subtitle} >  <b>Peptide sequence: </b> 9-26AA long peptide sequences (one per line).</div>
        </div>
          <br />
          <button style={styles.button} onClick={() => postData(input_allele,input_gene,input_seq, this)}>Analyze</button>
        <DataDisplay scores={this.state.scores}/>
      </div>);
  }
}

DataEntry.styles = {
  subtitle: {
    margin: "0px 1px 10px 1px",
    width: "92%"
  },
  pane: {
    display:"inline-block",
    margin: "20px auto 5px auto",
    width: "30%",
  },
  textarea: {
    width: "90%",
    height: "300px",
    margin: "0px 0px 5px 0px",
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
    return (<div>
      <div className="label" style={styles.subtitle} > <b>MARIA predicted presentation scores (assuming gene expression of B-cell lymphoma):</b></div>
      <textarea className="right_pane" style={pane_style} value={format_scores}></textarea>
      </div>);
  }
}

DataDisplay.styles = {
  subtitle: {
    margin: "10px 1px 3px 1px",
    width: "92%"
  },
  pane: {
    margin: "5px 1px 5px 1px",
    width: "70%",
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
        <h1 style={styles.h1}><b>M</b>HC <b>A</b>nalysis with <b>R</b>ecurrent <b>I</b>tergrated <b>A</b>rchitecture </h1>
        <h2 style={styles.h2}>A deep learning based tool for predicting MHC-II peptide presentation</h2>
        <DataEntry />
      </div>);
  }
}

App.styles = {
  body: {
    fontSize: "14px",
    fontFamily: "helvetica neue, helvetica, sans-serif",
    width: "1200px",
    margin: "20px auto"
  },
  h1: {
    fontSize: "2.0em",
    fontWeight: "400",
    margin: "3px 0px",
  },
  h2: {
    fontSize: "1.3em",
    fontWeight: "200",
    margin: "3px 0px"
  }
}

render(<App/>, document.getElementById('container'));
