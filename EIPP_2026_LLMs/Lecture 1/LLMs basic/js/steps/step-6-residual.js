// step-6-residual.js — Step 6: Residual Connections & Layer Norm
(function() {
  'use strict';
  var viz = window.TP.viz;
  var state = window.TP.state;

  var content = null;
  var flowContainer = null;
  var distContainer = null;
  var formulaContainer = null;

  function getData() {
    if (!state.result) return null;
    // Pick a representative token (last for DNA, mask position for protein)
    var tokenIdx = state.mode === 'dna'
      ? state.sequence.length - 1
      : state.result.outputIdx;
    return {
      seq: state.sequence,
      mode: state.mode,
      tokenIdx: tokenIdx,
      // Inputs to residual: positioned embeddings (X) and multiHeadOutput
      positionedEmbeddings: state.result.positionedEmbeddings,
      multiHeadOutput: state.result.multiHeadOutput,
      postResidual1: state.result.postResidual1,
      postLayerNorm1: state.result.postLayerNorm1,
      // Second residual: postLayerNorm1 + ffnOutput
      ffnOutput: state.result.ffnOutput,
      postResidual2: state.result.postResidual2,
      postLayerNorm2: state.result.postLayerNorm2
    };
  }

  function renderFlow() {
    if (!flowContainer) return;
    flowContainer.innerHTML = '';

    var card = document.createElement('div');
    card.className = 'tp-card';

    var title = document.createElement('div');
    title.style.cssText = 'font-family: var(--font-display); font-size: 0.9rem; font-weight: 600; margin-bottom: var(--space-sm);';
    title.textContent = 'Residual + LayerNorm Flow (Sub-layer 1: Attention)';
    card.appendChild(title);

    var hint = document.createElement('p');
    hint.style.cssText = 'font-size: 0.82rem; color: var(--text-muted); margin-bottom: var(--space-sm);';
    hint.textContent = 'Press Play to see how the residual connection preserves the original signal while adding the attention output, then LayerNorm stabilizes the result.';
    card.appendChild(hint);

    // Flow diagram container
    var diagramArea = document.createElement('div');
    diagramArea.style.cssText = 'display: flex; flex-direction: column; align-items: center; gap: 6px; padding: var(--space-md) 0; opacity: 0; transition: opacity 0.4s ease;';

    function makeBox(label, cls) {
      var box = document.createElement('div');
      box.style.cssText = 'padding: 8px 20px; border-radius: 6px; font-family: var(--font-mono); font-size: 0.78rem; font-weight: 600; text-align: center; transition: all 0.4s ease;';
      if (cls === 'input') {
        box.style.background = '#eef6f8';
        box.style.border = '1px solid #c8dfe6';
        box.style.color = '#3a6d7d';
      } else if (cls === 'sublayer') {
        box.style.background = 'var(--accent-light)';
        box.style.border = '1px solid var(--accent)';
        box.style.color = 'var(--accent)';
      } else if (cls === 'add') {
        box.style.background = '#f0f0e8';
        box.style.border = '1px solid #bbb';
        box.style.color = '#555';
        box.style.borderRadius = '50%';
        box.style.width = '36px';
        box.style.height = '36px';
        box.style.padding = '0';
        box.style.display = 'flex';
        box.style.alignItems = 'center';
        box.style.justifyContent = 'center';
        box.style.fontSize = '1.2rem';
      } else if (cls === 'norm') {
        box.style.background = '#e8f5e9';
        box.style.border = '1px solid #81c784';
        box.style.color = '#2e7d32';
      } else if (cls === 'output') {
        box.style.background = '#fff3e0';
        box.style.border = '1px solid #ffb74d';
        box.style.color = '#e65100';
      }
      box.textContent = label;
      return box;
    }

    function makeArrow() {
      var arrow = document.createElement('div');
      arrow.style.cssText = 'font-size: 1rem; color: var(--text-muted); line-height: 1;';
      arrow.textContent = '\u2193';
      return arrow;
    }

    // Build flow: Input X → fork → [skip, attention] → + → LayerNorm → Output
    var inputBox = makeBox('X (positioned embeddings)', 'input');
    var forkLabel = document.createElement('div');
    forkLabel.style.cssText = 'display: flex; gap: 40px; align-items: flex-start;';

    var skipPath = document.createElement('div');
    skipPath.style.cssText = 'display: flex; flex-direction: column; align-items: center; gap: 6px;';
    var skipLabel = document.createElement('div');
    skipLabel.style.cssText = 'font-size: 0.72rem; color: var(--text-muted); font-style: italic;';
    skipLabel.textContent = 'skip (identity)';
    var skipArrow = document.createElement('div');
    skipArrow.style.cssText = 'width: 2px; height: 40px; background: #bbb; border-radius: 1px;';
    skipPath.appendChild(skipLabel);
    skipPath.appendChild(skipArrow);

    var attnPath = document.createElement('div');
    attnPath.style.cssText = 'display: flex; flex-direction: column; align-items: center; gap: 6px;';
    var attnBox = makeBox('Multi-Head Attention', 'sublayer');
    var attnArrow = makeArrow();
    attnPath.appendChild(attnBox);
    attnPath.appendChild(attnArrow);

    forkLabel.appendChild(skipPath);
    forkLabel.appendChild(attnPath);

    var addBox = makeBox('+', 'add');
    var normBox = makeBox('LayerNorm', 'norm');
    var outputBox = makeBox('x\u2032 (normalized)', 'output');

    diagramArea.appendChild(inputBox);
    diagramArea.appendChild(makeArrow());
    diagramArea.appendChild(forkLabel);
    diagramArea.appendChild(addBox);
    diagramArea.appendChild(makeArrow());
    diagramArea.appendChild(normBox);
    diagramArea.appendChild(makeArrow());
    diagramArea.appendChild(outputBox);

    // Play animation
    viz.createPlayButton(card, function(done) {
      diagramArea.style.opacity = '1';
      // Phase 1: highlight input
      inputBox.style.boxShadow = '0 0 8px var(--accent)';
      setTimeout(function() {
        // Phase 2: highlight fork paths
        skipArrow.style.background = 'var(--accent)';
        attnBox.style.boxShadow = '0 0 8px var(--accent)';
      }, 600);
      setTimeout(function() {
        // Phase 3: addition
        addBox.style.boxShadow = '0 0 8px #888';
        addBox.style.background = '#e0e0d8';
      }, 1200);
      setTimeout(function() {
        // Phase 4: normalization
        normBox.style.boxShadow = '0 0 8px #4caf50';
      }, 1800);
      setTimeout(function() {
        // Phase 5: output
        outputBox.style.boxShadow = '0 0 8px #ff9800';
        done();
      }, 2400);
    });

    card.appendChild(diagramArea);
    flowContainer.appendChild(card);
  }

  function renderDistributions() {
    if (!distContainer) return;
    distContainer.innerHTML = '';
    var data = getData();
    if (!data) return;

    var title = document.createElement('div');
    title.className = 'tp-viz-title';
    title.textContent = 'Before vs After LayerNorm (token ' + data.tokenIdx + ': ' + data.seq[data.tokenIdx] + ')';
    distContainer.appendChild(title);

    var hint = document.createElement('p');
    hint.style.cssText = 'font-size: 0.82rem; color: var(--text-muted); margin-bottom: var(--space-sm);';
    hint.textContent = 'Before: the residual sum (X + Attention) has variable mean and variance. After LayerNorm: mean \u2248 0, variance \u2248 1.';
    distContainer.appendChild(hint);

    // Sub-layer 1: attention residual
    var before1 = data.postResidual1[data.tokenIdx];
    var after1 = data.postLayerNorm1[data.tokenIdx];

    var labels = [];
    for (var i = 0; i < TP.D_MODEL; i++) labels.push('d' + i);

    var dist1 = document.createElement('div');
    dist1.style.marginBottom = 'var(--space-md)';
    viz.renderDistributionBars(before1, after1, dist1, { labels: labels });
    distContainer.appendChild(dist1);

    // Sub-layer 2: FFN residual
    var subTitle2 = document.createElement('div');
    subTitle2.className = 'tp-viz-title';
    subTitle2.style.marginTop = 'var(--space-md)';
    subTitle2.textContent = 'Sub-layer 2: FFN Residual + LayerNorm';
    distContainer.appendChild(subTitle2);

    var before2 = data.postResidual2[data.tokenIdx];
    var after2 = data.postLayerNorm2[data.tokenIdx];

    var dist2 = document.createElement('div');
    viz.renderDistributionBars(before2, after2, dist2, { labels: labels });
    distContainer.appendChild(dist2);
  }

  function renderFormula() {
    if (!formulaContainer) return;
    formulaContainer.innerHTML = '';

    viz.renderFormulaBlock([
      'Residual: x\u2032 = x + SubLayer(x)',
      'LayerNorm(x) = \u03B3 \u2299 ((x \u2212 \u03BC) / \u221A(\u03C3\u00B2 + \u03B5)) + \u03B2',
      '\u03BC = mean(x), \u03C3\u00B2 = var(x)  across d_model=' + TP.D_MODEL + ' dimensions',
      '\u03B3, \u03B2 \u2208 \u211D<sup>d_model</sup> (learned scale & shift), \u03B5 = 10<sup>\u22125</sup>',
      'Full block: x\u2032 = LN(x + MHA(x)), x\u2033 = LN(x\u2032 + FFN(x\u2032))'
    ], formulaContainer);
  }

  function buildContent() {
    content = document.getElementById('step-6-content');
    content.innerHTML = '';

    // --- Text ---
    var p1 = document.createElement('p');
    p1.textContent = 'Two critical components make deep transformers trainable. Without them, stacking many layers would lead to vanishing or exploding gradients and training instability.';
    content.appendChild(p1);

    var p2 = document.createElement('p');
    p2.textContent = 'Residual connections (skip connections) add the input of each sub-layer directly to its output: x\u2032 = x + SubLayer(x). This ensures information can flow directly through the network \u2014 each layer only needs to learn a small correction ("delta"), which is much easier to optimize.';
    content.appendChild(p2);

    var p3 = document.createElement('p');
    p3.textContent = 'Layer normalization normalizes activations across the feature dimension (not the batch dimension), keeping values in a well-behaved range. The full transformer block applies two residual + LayerNorm pairs: one around attention, one around the FFN.';
    content.appendChild(p3);

    // --- Biology analogy ---
    viz.renderAnalogyBox(
      'Residual connections are like homeostatic feedback mechanisms in cells. No matter how many signaling cascades a signal passes through, the cell maintains a baseline state. If a particular transformation "malfunctions" (produces near-zero gradients), the residual ensures the original signal persists \u2014 just as homeostasis preserves baseline cellular function despite perturbations.',
      content
    );

    // --- Animated flow diagram ---
    flowContainer = document.createElement('div');
    flowContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(flowContainer);

    // --- Distribution bars ---
    distContainer = document.createElement('div');
    distContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(distContainer);

    // --- Formula block ---
    formulaContainer = document.createElement('div');
    formulaContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(formulaContainer);

    // Render all
    renderFlow();
    renderDistributions();
    renderFormula();
  }

  window.TP.steps[6] = {
    init: function() {},
    render: function() {
      buildContent();
    }
  };
})();
