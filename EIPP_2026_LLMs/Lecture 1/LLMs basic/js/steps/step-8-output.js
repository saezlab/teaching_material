// step-8-output.js — Step 8: Output & Prediction
(function() {
  'use strict';
  var viz = window.TP.viz;
  var state = window.TP.state;

  var content = null;
  var barContainer = null;
  var resultContainer = null;
  var captionContainer = null;
  var formulaContainer = null;

  // Amino acid physicochemical groups for caption
  var AA_GROUPS = {
    hydrophobic: 'LIVAFW M',
    charged_pos: 'KR',
    charged_neg: 'DE',
    polar: 'STNQ',
    special: 'GPHY C'
  };

  function getAAGroup(aa) {
    if ('LIVAFWM'.indexOf(aa) >= 0) return 'hydrophobic';
    if ('KR'.indexOf(aa) >= 0) return 'positively charged';
    if ('DE'.indexOf(aa) >= 0) return 'negatively charged';
    if ('STNQ'.indexOf(aa) >= 0) return 'polar';
    if (aa === 'G') return 'small/flexible (glycine)';
    if (aa === 'P') return 'rigid (proline)';
    if (aa === 'H') return 'polar/charged (histidine)';
    if (aa === 'Y') return 'aromatic/polar (tyrosine)';
    if (aa === 'C') return 'thiol-containing (cysteine)';
    return 'standard';
  }

  // DNA complement map
  var DNA_COMPLEMENT = { A: 'T', T: 'A', G: 'C', C: 'G' };

  function getData() {
    if (!state.result) return null;
    return {
      seq: state.sequence,
      mode: state.mode,
      logits: state.result.logits,
      probabilities: state.result.probabilities,
      predictedToken: state.result.predictedToken,
      outputIdx: state.result.outputIdx,
      tokenIds: state.result.tokenIds
    };
  }

  function getVocabLabels() {
    return state.mode === 'dna' ? TP.DNA_VOCAB_REV.slice() : TP.PROTEIN_VOCAB_REV.slice();
  }

  function renderBarChart() {
    if (!barContainer) return;
    barContainer.innerHTML = '';
    var data = getData();
    if (!data) return;

    var labels = getVocabLabels();
    var topIdx = TP.argmax(data.probabilities);

    viz.renderBarChart(labels, data.probabilities, barContainer, {
      highlightIdx: topIdx,
      maxHeight: state.mode === 'dna' ? 180 : 140,
      animate: true
    });
  }

  function renderResult() {
    if (!resultContainer) return;
    resultContainer.innerHTML = '';
    var data = getData();
    if (!data) return;

    var tokens = data.seq.split('');
    var topIdx = TP.argmax(data.probabilities);
    var topProb = data.probabilities[topIdx];

    // Build the result token row
    var title = document.createElement('div');
    title.className = 'tp-viz-title';

    if (data.mode === 'dna') {
      title.textContent = 'Prediction: next token after "' + data.seq + '"';
    } else {
      title.textContent = 'Prediction: fill [MASK] at position ' + data.outputIdx;
    }
    resultContainer.appendChild(title);

    // Build display tokens
    var displayTokens, displayIds;
    if (data.mode === 'dna') {
      // Append predicted token
      displayTokens = tokens.concat([data.predictedToken]);
      displayIds = data.tokenIds.concat([topIdx]);
    } else {
      // Replace mask position with prediction
      displayTokens = tokens.slice();
      displayIds = data.tokenIds.slice();
      displayTokens[data.outputIdx] = data.predictedToken;
      displayIds[data.outputIdx] = topIdx;
    }

    var rowWrapper = document.createElement('div');

    viz.renderTokenRow(displayTokens, rowWrapper, {
      ids: displayIds,
      mode: data.mode
    });

    // Add highlight on the predicted position
    var tokenBlocks = rowWrapper.querySelectorAll('.tp-token-block');
    var highlightIdx = data.mode === 'dna' ? displayTokens.length - 1 : data.outputIdx;
    if (tokenBlocks[highlightIdx]) {
      tokenBlocks[highlightIdx].style.outline = '3px solid var(--accent)';
      tokenBlocks[highlightIdx].style.outlineOffset = '2px';
    }

    resultContainer.appendChild(rowWrapper);

    // Probability badge
    var badge = document.createElement('div');
    badge.style.cssText = 'display: inline-block; margin-top: 8px; padding: 4px 12px; border-radius: 12px; background: var(--accent-light); color: var(--accent); font-family: var(--font-mono); font-size: 0.8rem; font-weight: 600;';
    badge.textContent = data.predictedToken + ': ' + (topProb * 100).toFixed(1) + '% confidence';
    resultContainer.appendChild(badge);
  }

  function renderCaption() {
    if (!captionContainer) return;
    captionContainer.innerHTML = '';
    var data = getData();
    if (!data) return;

    var topIdx = TP.argmax(data.probabilities);
    var topProb = data.probabilities[topIdx];
    var predicted = data.predictedToken;

    var caption = document.createElement('div');
    caption.style.cssText = 'padding: var(--space-md); background: var(--accent-light); border-radius: var(--border-radius); border-left: 3px solid var(--accent); margin-top: var(--space-md);';

    var text = '';
    if (data.mode === 'dna') {
      var lastBase = data.seq[data.seq.length - 1];
      var complement = DNA_COMPLEMENT[lastBase] || '?';
      if (predicted === complement) {
        text = 'The model predicts ' + predicted + ' with ' + (topProb * 100).toFixed(1) + '% confidence. This makes biological sense: ' + predicted + ' is the Watson-Crick complement of the last base ' + lastBase + '. The attention mechanism learned that in DNA, each base is most strongly associated with its complement (A\u2194T, G\u2194C), mimicking base-pairing rules.';
      } else {
        text = 'The model predicts ' + predicted + ' with ' + (topProb * 100).toFixed(1) + '% confidence. The complement of the last base ' + lastBase + ' would be ' + complement + '. The prediction reflects the learned patterns in the rigged weights, which capture complementary base-pairing relationships.';
      }
    } else {
      var maskedAA = data.seq[data.outputIdx];
      var maskedGroup = getAAGroup(maskedAA);
      var predictedGroup = getAAGroup(predicted);

      if (maskedGroup === predictedGroup || (maskedGroup.indexOf('hydrophobic') >= 0 && predictedGroup.indexOf('hydrophobic') >= 0)) {
        text = 'The model predicts ' + predicted + ' (' + predictedGroup + ') with ' + (topProb * 100).toFixed(1) + '% confidence to replace the masked ' + maskedAA + ' (' + maskedGroup + '). This is physicochemically sensible: the model predicts a residue from the same class, reflecting learned substitution patterns similar to BLOSUM62.';
      } else {
        text = 'The model predicts ' + predicted + ' (' + predictedGroup + ') with ' + (topProb * 100).toFixed(1) + '% confidence to replace the masked ' + maskedAA + ' (' + maskedGroup + '). While not from the same physicochemical class, the prediction reflects the surrounding sequence context captured by bidirectional attention.';
      }
    }

    var textEl = document.createElement('p');
    textEl.style.cssText = 'margin: 0; font-size: 0.88rem; line-height: 1.5;';
    textEl.textContent = text;
    caption.appendChild(textEl);
    captionContainer.appendChild(caption);
  }

  function renderFormula() {
    if (!formulaContainer) return;
    formulaContainer.innerHTML = '';

    var vocabSize = state.mode === 'dna' ? 4 : 20;

    viz.renderFormulaBlock([
      'logits = <b>h</b> \u00B7 W<sub>vocab</sub> + b<sub>vocab</sub> \u2208 \u211D<sup>|V|</sup>  (|V| = ' + vocabSize + ')',
      'P(token<sub>i</sub>) = softmax(logits)<sub>i</sub> = exp(logit<sub>i</sub>) / \u03A3<sub>j</sub> exp(logit<sub>j</sub>)',
      'Loss = \u2212log P(correct token)  (cross-entropy)',
      'Temperature: P(token<sub>i</sub>) = softmax(logits / T)<sub>i</sub>  &mdash;  T &lt; 1 \u2192 sharper, T &gt; 1 \u2192 more uniform'
    ], formulaContainer);
  }

  function buildContent() {
    content = document.getElementById('step-8-content');
    content.innerHTML = '';

    // --- Text ---
    var p1 = document.createElement('p');
    p1.textContent = 'The final transformer representation is projected into vocabulary space using a linear layer, producing one logit per vocabulary item. Softmax converts these logits into a probability distribution over all possible next tokens.';
    content.appendChild(p1);

    var p2 = document.createElement('p');
    if (state.mode === 'dna') {
      p2.textContent = 'In DNA (decoder) mode, we take the last token\u2019s representation and predict the next nucleotide. During training, this is done for every position: position i predicts position i+1. The training loss is cross-entropy: \u2212log P(correct token).';
    } else {
      p2.textContent = 'In protein (encoder) mode, we take the masked position\u2019s representation and predict the original residue. This is the masked language model (MLM) objective used by BERT and ESM. The training loss is cross-entropy: \u2212log P(correct token).';
    }
    content.appendChild(p2);

    // --- Biology analogy ---
    viz.renderAnalogyBox(
      'This is the functional output: the gene product, the final phenotype. All the prior computation \u2014 embeddings as amino acid properties, attention as residue interactions, FFN as post-translational processing \u2014 culminates in a prediction. In DNA mode, it\u2019s like predicting the next base in a replication fork given the template strand. In protein mode, it\u2019s like predicting which residue would be most structurally and chemically compatible at a gap position.',
      content
    );

    // --- Bar chart ---
    var barTitle = document.createElement('h3');
    barTitle.textContent = 'Prediction Probabilities';
    barTitle.style.cssText = 'margin-top: var(--space-xl); margin-bottom: var(--space-sm);';
    content.appendChild(barTitle);

    var barHint = document.createElement('p');
    barHint.style.cssText = 'font-size: 0.85rem; color: var(--text-muted); margin-bottom: var(--space-sm);';
    barHint.textContent = state.mode === 'dna'
      ? 'Probability for each nucleotide as the next token. The highlighted bar is the top prediction.'
      : 'Probability for each amino acid at the masked position. The highlighted bar is the top prediction.';
    content.appendChild(barHint);

    barContainer = document.createElement('div');
    content.appendChild(barContainer);

    // --- Result token row ---
    resultContainer = document.createElement('div');
    resultContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(resultContainer);

    // --- Biological sensibility caption ---
    captionContainer = document.createElement('div');
    content.appendChild(captionContainer);

    // --- Formula block ---
    formulaContainer = document.createElement('div');
    formulaContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(formulaContainer);

    // Render all
    renderBarChart();
    renderResult();
    renderCaption();
    renderFormula();
  }

  window.TP.steps[8] = {
    init: function() {},
    render: function() {
      buildContent();
    }
  };
})();
