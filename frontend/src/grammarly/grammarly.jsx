import React from 'react';
import ReactDOM from 'react-dom';
import { List } from 'react-virtualized';

import { checkLipinski, checkProps, checkValidity, getImage } from './service';

import './grammarly.css';

function submit(smiles) {
  checkValidity(smiles).then(function (result) {
    console.log(result);
  });
}
window.submit = submit;

function props(smiles) {
  checkProps(smiles).then(function (result) {
    console.log(result);
  });
}
window.props = props;

function lipinski(smiles) {
  checkLipinski(smiles).then(function (result) {
    console.log(result);
  });
}
window.lipinski = lipinski;

function gImage(smiles) {
  getImage(smiles).then(function (result) {
    console.log(result);
  });
}
window.gImage = gImage;

function Grammarly() {
  return (
    <div className="main_content">
      {/* logo */}
      <div className="logo">
        <img src="images/smiles.png" title="SMILES logo" width="50" height="50"></img>
        <h1>Grammarly4Smiles</h1>
      </div>

      {/* first chatbox */}
      <List rowRenderer={Row} rowCount={3} rowHeight={160} height="300" width="300" />
      <form>
        <label>
          <input type="text" name="information" placeholder="Type here ..." />
        </label>
        <input type="submit"></input>
      </form>
    </div>
  );
}

function Row({ index }) {
  return [
    <h4>Smiling Assistant | Customer Service</h4>,
    <div className="Assistant1">
      <img src="images/smiles.png" title="SMILES logo" width="45" height="45"></img>
      <h3>Hi! Smiling here! Please enter your beautiful SMILES!</h3>
    </div>,
    <p className="time">10:52 AM</p>,
  ];
}

export default Grammarly;

// scrolling task and how to hook up a data structure pass the messages and whether it is the th
