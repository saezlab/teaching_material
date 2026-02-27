// step-4-attention.js — Step 4: Self-Attention Deep Dive
(function() {
  'use strict';
  var viz = window.TP.viz;
  var state = window.TP.state;

  var content = null;
  var stageContainer = null;
  var tokenRowContainer = null;
  var formulaContainer = null;
  var tokenRowEl = null; // the actual .tp-token-row element for arrows

  function getData() {
    if (!state.result || !state.result.qkvPerHead) return null;
    var head = state.result.qkvPerHead[0]; // Head 0 only
    return {
      seq: state.sequence,
      mode: state.mode,
      tokens: state.sequence.split(''),
      tokenIds: state.result.tokenIds,
      Q: head.Q,
      K: head.K,
      V: head.V,
      scores: head.scores,
      scaledScores: head.scaledScores,
      mask: head.mask,
      attention: head.attention
    };
  }

  var getTokenLabels = TP.getTokenLabels;

  function getDkLabels() {
    return ['k0', 'k1'];
  }

  // Build the masked/scaled scores matrix for display
  // In DNA mode, masked cells show as NaN (we'll render them specially)
  function getMaskedScaledMatrix(data) {
    var n = data.tokens.length;
    var mat = [];
    for (var i = 0; i < n; i++) {
      var row = [];
      for (var j = 0; j < n; j++) {
        if (data.mask[i][j] === 1) {
          row.push(null); // masked
        } else {
          row.push(data.scaledScores[i][j]);
        }
      }
      mat.push(row);
    }
    return mat;
  }

  // Render a matrix that may have null (masked) cells
  function renderMaskedMatrix(data, matrix, container, options) {
    options = options || {};
    container.innerHTML = '';

    var rows = matrix.length;
    var cols = matrix[0].length;
    var tokenLabels = getTokenLabels(data.tokens);

    // Find min/max of non-null values
    var min = 0, max = 0;
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < cols; j++) {
        if (matrix[i][j] !== null) {
          if (matrix[i][j] < min) min = matrix[i][j];
          if (matrix[i][j] > max) max = matrix[i][j];
        }
      }
    }

    if (options.title) {
      var titleEl = document.createElement('div');
      titleEl.className = 'tp-viz-title';
      titleEl.textContent = options.title;
      container.appendChild(titleEl);
    }

    var wrapper = document.createElement('div');
    wrapper.className = 'tp-matrix-container';

    var grid = document.createElement('div');
    grid.className = 'tp-matrix';
    grid.style.gridTemplateColumns = 'repeat(' + (cols + 1) + ', auto)';

    // Column labels
    var corner = document.createElement('div');
    corner.className = 'tp-matrix-label';
    grid.appendChild(corner);
    for (var j = 0; j < cols; j++) {
      var cl = document.createElement('div');
      cl.className = 'tp-matrix-label';
      cl.textContent = tokenLabels[j];
      grid.appendChild(cl);
    }

    // Data rows
    for (var i = 0; i < rows; i++) {
      var rl = document.createElement('div');
      rl.className = 'tp-matrix-label';
      rl.textContent = tokenLabels[i];
      grid.appendChild(rl);

      for (var j = 0; j < cols; j++) {
        var cell = document.createElement('div');
        cell.className = 'tp-matrix-cell';

        if (matrix[i][j] === null) {
          // Masked cell
          cell.style.backgroundColor = '#2a2a2a';
          cell.style.color = '#666';
          cell.textContent = '\u2013';
        } else {
          var v = matrix[i][j];
          if (options.scale === 'greenScale') {
            cell.style.backgroundColor = viz.valueToGreenScale(v);
            cell.style.color = v > 0.55 ? '#fff' : '#333';
          } else {
            var range = Math.max(Math.abs(min), Math.abs(max));
            if (range === 0) range = 1;
            cell.style.backgroundColor = viz.valueToBlueWhiteRed(v, min, max);
            cell.style.color = Math.abs(v / range) > 0.6 ? '#fff' : '#333';
          }
          cell.textContent = v.toFixed(2);
        }

        grid.appendChild(cell);
      }
    }

    wrapper.appendChild(grid);
    container.appendChild(wrapper);
  }

  // Render the attention heatmap with clickable cells
  function renderAttentionHeatmap(data, container) {
    container.innerHTML = '';

    var tokenLabels = getTokenLabels(data.tokens);
    var attn = data.attention;

    viz.renderMatrix(attn, container, {
      title: 'Stage 4: Attention Weights (softmax) \u2014 click a row to see arrows',
      rowLabels: tokenLabels,
      colLabels: tokenLabels,
      scale: 'greenScale',
      onClick: function(row, col, val) {
        // Draw arrows for this row (source token = row)
        if (tokenRowEl) {
          viz.renderArrows(tokenRowEl, attn[row], row);
        }
      }
    });
  }

  function renderTokenRow() {
    if (!tokenRowContainer) return;
    tokenRowContainer.innerHTML = '';
    var data = getData();
    if (!data) return;

    // Add some spacing and a label above the token row
    var label = document.createElement('div');
    label.style.cssText = 'font-size: 0.82rem; color: var(--text-muted); margin-bottom: var(--space-xs); margin-top: var(--space-sm);';
    label.textContent = 'Click any cell in the attention heatmap to see arrows:';
    tokenRowContainer.appendChild(label);

    // Wrapper needs position:relative for SVG arrows
    var arrowWrapper = document.createElement('div');
    arrowWrapper.style.cssText = 'position: relative; padding-top: 90px;'; // space for arrows above

    var maskIdx = (data.mode === 'protein') ? Math.floor(data.seq.length / 2) : undefined;
    tokenRowEl = viz.renderTokenRow(data.tokens, arrowWrapper, {
      ids: data.tokenIds,
      maskIdx: maskIdx,
      mode: data.mode
    });

    tokenRowContainer.appendChild(arrowWrapper);
  }

  function renderFormula() {
    if (!formulaContainer) return;
    formulaContainer.innerHTML = '';

    viz.renderFormulaBlock([
      'Q = XW<sub>Q</sub>, K = XW<sub>K</sub>, V = XW<sub>V</sub>  (W<sub>Q</sub>, W<sub>K</sub>, W<sub>V</sub> \u2208 \u211D<sup>d_model \u00D7 d_k</sup>)',
      'scores = QK<sup>T</sup> \u2208 \u211D<sup>n \u00D7 n</sup>',
      'scaled = scores / \u221Ad<sub>k</sub>  (scaling prevents softmax saturation)',
      (state.mode === 'dna' ? 'masked: scaled + M, where M<sub>ij</sub> = \u2212\u221E if j > i (causal mask)' : 'No mask in encoder/BERT mode \u2014 all positions attend to all others'),
      'attn = softmax(scaled) per row  (each row sums to 1)',
      'output = attn \u00B7 V \u2208 \u211D<sup>n \u00D7 d_k</sup>'
    ], formulaContainer);
  }

  function runAnimation(data) {
    if (!stageContainer) return;
    stageContainer.innerHTML = '';

    var tokenLabels = getTokenLabels(data.tokens);
    var dkLabels = getDkLabels();

    // Create 4 stage divs, initially hidden
    var stages = [];
    for (var s = 0; s < 4; s++) {
      var stage = document.createElement('div');
      stage.style.cssText = 'opacity: 0; transition: opacity 0.4s ease; margin-bottom: var(--space-md);';
      stageContainer.appendChild(stage);
      stages.push(stage);
    }

    // Stage 1: Q, K, V matrices side-by-side
    var s1 = stages[0];
    var s1Title = document.createElement('div');
    s1Title.className = 'tp-viz-title';
    s1Title.textContent = 'Stage 1: Query, Key, Value Projections (Head 0)';
    s1.appendChild(s1Title);

    var s1Row = document.createElement('div');
    s1Row.style.cssText = 'display: flex; gap: var(--space-md); flex-wrap: wrap;';

    var qDiv = document.createElement('div');
    viz.renderMatrix(data.Q, qDiv, {
      title: 'Q (Query)',
      rowLabels: tokenLabels,
      colLabels: dkLabels,
      scale: 'blueWhiteRed'
    });
    s1Row.appendChild(qDiv);

    var kDiv = document.createElement('div');
    viz.renderMatrix(data.K, kDiv, {
      title: 'K (Key)',
      rowLabels: tokenLabels,
      colLabels: dkLabels,
      scale: 'blueWhiteRed'
    });
    s1Row.appendChild(kDiv);

    var vDiv = document.createElement('div');
    viz.renderMatrix(data.V, vDiv, {
      title: 'V (Value)',
      rowLabels: tokenLabels,
      colLabels: dkLabels,
      scale: 'blueWhiteRed'
    });
    s1Row.appendChild(vDiv);

    s1.appendChild(s1Row);

    // Stage 2: QK^T raw scores
    var s2 = stages[1];
    viz.renderMatrix(data.scores, s2, {
      title: 'Stage 2: Raw Scores QK\u1D40 (' + data.tokens.length + '\u00D7' + data.tokens.length + ')',
      rowLabels: tokenLabels,
      colLabels: tokenLabels,
      scale: 'blueWhiteRed'
    });

    // Stage 3: Scaled + masked scores
    var s3 = stages[2];
    var maskedMatrix = getMaskedScaledMatrix(data);
    renderMaskedMatrix(data, maskedMatrix, s3, {
      title: 'Stage 3: Scaled (\u00F7\u221A' + TP.D_K + ') ' + (data.mode === 'dna' ? '+ Causal Mask' : '(no mask)'),
      scale: 'blueWhiteRed'
    });

    // Stage 4: Softmax attention heatmap (clickable)
    var s4 = stages[3];
    renderAttentionHeatmap(data, s4);

    // Animate stages sequentially
    var delay = 300;
    for (var i = 0; i < 4; i++) {
      (function(idx) {
        setTimeout(function() {
          stages[idx].style.opacity = '1';
        }, delay + idx * 500);
      })(i);
    }

    return stages;
  }

  function buildContent() {
    content = document.getElementById('step-4-content');
    content.innerHTML = '';

    // --- Text ---
    var p1 = document.createElement('p');
    p1.textContent = 'Self-attention is the core of the transformer. It allows each token to dynamically decide which other tokens in the sequence are relevant to it. Three projections are computed per token: Query (\u201CWhat am I looking for?\u201D), Key (\u201CWhat do I contain?\u201D), and Value (\u201CWhat information should I pass along?\u201D).';
    content.appendChild(p1);

    var p2 = document.createElement('p');
    p2.textContent = 'Attention works like a soft dictionary lookup: the Query of token i is compared against all Keys via dot product to produce a relevance score. These scores are normalized via softmax to sum to 1, then used as weights to compute a weighted average of Values. The scaling factor 1/\u221Ad_k prevents dot products from growing too large, which would push softmax into saturation.';
    content.appendChild(p2);

    var p3 = document.createElement('p');
    if (state.mode === 'dna') {
      p3.textContent = 'In decoder (DNA) mode, a causal mask sets scores for future positions to \u2212\u221E before softmax, ensuring the model can only attend to positions \u2264 the current position. This is shown below as blacked-out cells in the upper triangle.';
    } else {
      p3.textContent = 'In encoder (protein) mode, all positions can attend to all others \u2014 there is no causal mask. This bidirectional attention allows each residue to gather context from the entire sequence.';
    }
    content.appendChild(p3);

    // --- Biology analogy ---
    viz.renderAnalogyBox(
      'The attention weight matrix has the same shape as a residue contact map in structural biology: an n\u00D7n matrix where entry (i,j) measures how much position i \u201Cinteracts with\u201D position j. In direct coupling analysis (DCA/EVcouplings), you compute pairwise statistical couplings from MSA data. Self-attention learns analogous pairwise relationships, but parameterized and differentiable, trained end-to-end.',
      content
    );

    // --- Play button + stage container ---
    var playCard = document.createElement('div');
    playCard.className = 'tp-card';
    playCard.style.marginTop = 'var(--space-lg)';

    var playTitle = document.createElement('div');
    playTitle.style.cssText = 'font-family: var(--font-display); font-size: 0.9rem; font-weight: 600; margin-bottom: var(--space-sm);';
    playTitle.textContent = 'Attention Computation (Head 0)';
    playCard.appendChild(playTitle);

    var playHint = document.createElement('p');
    playHint.style.cssText = 'font-size: 0.82rem; color: var(--text-muted); margin-bottom: var(--space-sm);';
    playHint.textContent = 'Press Play to see the 4 stages of self-attention unfold. After Stage 4, click any cell in the attention heatmap to see arrows between tokens.';
    playCard.appendChild(playHint);

    viz.createPlayButton(playCard, function(done) {
      var data = getData();
      if (!data) { done(); return; }

      runAnimation(data);
      renderTokenRow();

      // Done after all 4 stages have appeared
      setTimeout(function() { done(); }, 2500);
    });

    stageContainer = document.createElement('div');
    stageContainer.style.marginTop = 'var(--space-md)';
    playCard.appendChild(stageContainer);

    content.appendChild(playCard);

    // --- Token row for arrows ---
    tokenRowContainer = document.createElement('div');
    content.appendChild(tokenRowContainer);

    // --- Formula block ---
    formulaContainer = document.createElement('div');
    formulaContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(formulaContainer);

    renderFormula();

    // Auto-run animation on first render so content isn't empty
    var data = getData();
    if (data) {
      runAnimation(data);
      renderTokenRow();
    }
  }

  window.TP.steps[4] = {
    init: function() {},
    render: function() {
      buildContent();
    }
  };
})();
