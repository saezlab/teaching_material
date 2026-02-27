// step-7-ffn.js — Step 7: Feed-Forward Network
(function() {
  'use strict';
  var viz = window.TP.viz;
  var state = window.TP.state;

  var content = null;
  var nodeContainer = null;
  var activationContainer = null;
  var formulaContainer = null;

  function getData() {
    if (!state.result) return null;
    var tokenIdx = state.mode === 'dna'
      ? state.sequence.length - 1
      : state.result.outputIdx;
    return {
      seq: state.sequence,
      mode: state.mode,
      tokenIdx: tokenIdx,
      ffnHidden: state.result.ffnHidden,
      ffnOutput: state.result.ffnOutput,
      postLayerNorm1: state.result.postLayerNorm1
    };
  }

  function renderNodeDiagram() {
    if (!nodeContainer) return;
    nodeContainer.innerHTML = '';
    var data = getData();
    if (!data) return;

    var title = document.createElement('div');
    title.className = 'tp-viz-title';
    title.textContent = 'FFN Node Diagram (token ' + data.tokenIdx + ': ' + data.seq[data.tokenIdx] + ')';
    nodeContainer.appendChild(title);

    var hint = document.createElement('p');
    hint.style.cssText = 'font-size: 0.82rem; color: var(--text-muted); margin-bottom: var(--space-sm);';
    hint.textContent = 'Press Play to see the signal propagate through the FFN. Dark hidden nodes are zeroed by ReLU (negative pre-activation values clipped to 0).';
    nodeContainer.appendChild(hint);

    // Get actual hidden activations for this token to determine ReLU zeroing
    var hiddenActivations = data.ffnHidden[data.tokenIdx];

    var diagramDiv = document.createElement('div');
    viz.renderNodeDiagram(TP.D_MODEL, 12, TP.D_MODEL, diagramDiv, {
      hiddenTotal: TP.D_FF,
      activations: hiddenActivations
    });
    nodeContainer.appendChild(diagramDiv);
  }

  function renderActivationSummary() {
    if (!activationContainer) return;
    activationContainer.innerHTML = '';
    var data = getData();
    if (!data) return;

    var hiddenActivations = data.ffnHidden[data.tokenIdx];
    var zeroCount = 0;
    var maxAct = 0;
    for (var i = 0; i < hiddenActivations.length; i++) {
      if (hiddenActivations[i] === 0) zeroCount++;
      if (hiddenActivations[i] > maxAct) maxAct = hiddenActivations[i];
    }

    var summary = document.createElement('div');
    summary.style.cssText = 'padding: var(--space-sm) var(--space-md); background: var(--card-bg); border: 1px solid var(--card-border); border-radius: var(--border-radius); font-family: var(--font-mono); font-size: 0.78rem; color: var(--text-secondary);';

    var lines = [];
    lines.push('Hidden layer: ' + TP.D_FF + ' neurons total');
    lines.push('ReLU zeroed: ' + zeroCount + '/' + TP.D_FF + ' neurons (' + ((zeroCount / TP.D_FF * 100).toFixed(0)) + '% sparse)');
    lines.push('Max activation: ' + maxAct.toFixed(3));
    lines.push('Active neurons: ' + (TP.D_FF - zeroCount) + ' \u2014 these encode patterns relevant to this token');

    for (var i = 0; i < lines.length; i++) {
      var line = document.createElement('div');
      line.textContent = lines[i];
      if (i > 0) line.style.marginTop = '2px';
      summary.appendChild(line);
    }

    activationContainer.appendChild(summary);
  }

  function renderFormula() {
    if (!formulaContainer) return;
    formulaContainer.innerHTML = '';

    viz.renderFormulaBlock([
      'FFN(x) = ReLU(x\u00B7W<sub>1</sub> + b<sub>1</sub>) \u00B7 W<sub>2</sub> + b<sub>2</sub>',
      'W<sub>1</sub> \u2208 \u211D<sup>' + TP.D_MODEL + '\u00D7' + TP.D_FF + '</sup>, W<sub>2</sub> \u2208 \u211D<sup>' + TP.D_FF + '\u00D7' + TP.D_MODEL + '</sup>',
      'd_ff = 4 \u00B7 d_model = ' + TP.D_FF + '  (expansion factor = 4)',
      'FFN params: 2\u00B7d_model\u00B7d_ff + d_ff + d_model = ' + (2 * TP.D_MODEL * TP.D_FF + TP.D_FF + TP.D_MODEL),
      'Modern variant (GLU): FFN(x) = (xW<sub>1</sub> \u2299 \u03C3(xW<sub>gate</sub>))W<sub>2</sub>'
    ], formulaContainer);
  }

  function buildContent() {
    content = document.getElementById('step-7-content');
    content.innerHTML = '';

    // --- Text ---
    var p1 = document.createElement('p');
    p1.textContent = 'After attention has mixed information across positions, each token\u2019s representation passes independently through a two-layer fully connected network. "Independently" is key: unlike attention, the FFN applies the exact same transformation to every position with no cross-position interaction.';
    content.appendChild(p1);

    var p2 = document.createElement('p');
    p2.textContent = 'The FFN expands the representation from d_model=' + TP.D_MODEL + ' to d_ff=' + TP.D_FF + ' (4\u00D7 expansion), applies ReLU non-linearity, then compresses back to d_model=' + TP.D_MODEL + '. Research suggests FFNs act as "key-value memories": different neurons in the expanded layer activate for different input patterns, storing factual associations learned during training.';
    content.appendChild(p2);

    // --- Biology analogy ---
    viz.renderAnalogyBox(
      'Think of post-translational modifications. After the protein is assembled and has folded (attention gathered the right spatial context), individual residues get modified: phosphorylation, glycosylation, ubiquitination. These modifications are position-specific and context-dependent, but they happen to each residue independently \u2014 just like the FFN processes each token independently.',
      content
    );

    // --- Node diagram ---
    nodeContainer = document.createElement('div');
    nodeContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(nodeContainer);

    // --- Activation summary ---
    activationContainer = document.createElement('div');
    activationContainer.style.marginTop = 'var(--space-md)';
    content.appendChild(activationContainer);

    // --- Formula block ---
    formulaContainer = document.createElement('div');
    formulaContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(formulaContainer);

    // Render all
    renderNodeDiagram();
    renderActivationSummary();
    renderFormula();
  }

  window.TP.steps[7] = {
    init: function() {},
    render: function() {
      buildContent();
    }
  };
})();
