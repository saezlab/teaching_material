// step-3-positional.js — Step 3: Positional Encoding
(function() {
  'use strict';
  var viz = window.TP.viz;
  var state = window.TP.state;

  var content = null;
  var canvasContainer = null;
  var beforeAfterContainer = null;
  var formulaContainer = null;

  function getData() {
    if (!state.result) return null;
    return {
      seq: state.sequence,
      mode: state.mode,
      tokenIds: state.result.tokenIds,
      embeddings: state.result.embeddings,
      positionalEncodings: state.result.positionalEncodings,
      positionedEmbeddings: state.result.positionedEmbeddings
    };
  }

  var getDimLabels = TP.getDimLabels;

  // Find duplicate token positions (same token ID at different positions)
  function findDuplicatePositions(tokenIds) {
    var seen = {};
    var duplicates = [];
    for (var i = 0; i < tokenIds.length; i++) {
      var id = tokenIds[i];
      if (seen[id] !== undefined) {
        if (duplicates.indexOf(seen[id]) < 0) {
          duplicates.push(seen[id]);
        }
        duplicates.push(i);
      } else {
        seen[id] = i;
      }
    }
    return duplicates;
  }

  function renderCanvas() {
    if (!canvasContainer) return;
    var data = getData();
    if (!data || data.seq.length < 4) {
      canvasContainer.innerHTML = '';
      return;
    }

    viz.renderSinusoidalCanvas(data.seq.length, TP.D_MODEL, canvasContainer, {
      width: 600,
      height: 260
    });
  }

  function renderBeforeAfter() {
    if (!beforeAfterContainer) return;
    var data = getData();
    if (!data || data.seq.length < 4) {
      beforeAfterContainer.innerHTML = '';
      return;
    }

    beforeAfterContainer.innerHTML = '';

    var tokens = data.seq.split('');
    var duplicates = findDuplicatePositions(data.tokenIds);
    var dimLabels = getDimLabels();

    // Create side-by-side layout
    var wrapper = document.createElement('div');
    wrapper.style.cssText = 'display: flex; gap: var(--space-lg); flex-wrap: wrap; align-items: flex-start;';

    // Before matrix (raw embeddings)
    var beforeDiv = document.createElement('div');
    beforeDiv.style.cssText = 'flex: 1; min-width: 280px;';

    var beforeRowLabels = [];
    for (var i = 0; i < tokens.length; i++) {
      beforeRowLabels.push(i + ':' + tokens[i]);
    }

    viz.renderMatrix(data.embeddings, beforeDiv, {
      title: 'Before: Raw Embeddings',
      rowLabels: beforeRowLabels,
      colLabels: dimLabels,
      highlightRows: duplicates,
      scale: 'blueWhiteRed'
    });

    wrapper.appendChild(beforeDiv);

    // After matrix (embeddings + PE)
    var afterDiv = document.createElement('div');
    afterDiv.style.cssText = 'flex: 1; min-width: 280px;';

    var afterRowLabels = [];
    for (var i = 0; i < tokens.length; i++) {
      afterRowLabels.push(i + ':' + tokens[i]);
    }

    viz.renderMatrix(data.positionedEmbeddings, afterDiv, {
      title: 'After: Embeddings + PE',
      rowLabels: afterRowLabels,
      colLabels: dimLabels,
      highlightRows: duplicates,
      scale: 'blueWhiteRed'
    });

    wrapper.appendChild(afterDiv);
    beforeAfterContainer.appendChild(wrapper);

    // Annotation about identical tokens diverging
    if (duplicates.length > 0) {
      var note = document.createElement('div');
      note.style.cssText = 'font-size: 0.8rem; color: var(--accent); margin-top: 8px; font-style: italic;';

      // Group duplicates by token ID
      var groups = {};
      for (var i = 0; i < data.tokenIds.length; i++) {
        var id = data.tokenIds[i];
        if (!groups[id]) groups[id] = [];
        groups[id].push(i);
      }
      var parts = [];
      for (var id in groups) {
        if (groups[id].length > 1) {
          var token = tokens[groups[id][0]];
          parts.push('\u201C' + token + '\u201D (positions ' + groups[id].join(', ') + ')');
        }
      }
      note.textContent = 'Highlighted rows: ' + parts.join('; ') + ' \u2014 identical in the \u201CBefore\u201D matrix but now different after positional encoding is added.';
      beforeAfterContainer.appendChild(note);
    }
  }

  function renderFormula() {
    if (!formulaContainer) return;
    var data = getData();
    if (!data || data.seq.length < 4) return;

    viz.renderFormulaBlock([
      'PE(pos, 2i) = sin(pos / 10000<sup>2i/d_model</sup>)',
      'PE(pos, 2i+1) = cos(pos / 10000<sup>2i/d_model</sup>)',
      '<b>x</b> = <b>x</b><sub>embed</sub> + PE(pos)',
      'Wavelengths range from 2\u03C0 to 10000\u00B72\u03C0 across dimensions'
    ], formulaContainer);
  }

  function buildContent() {
    content = document.getElementById('step-3-content');
    content.innerHTML = '';

    // --- Text ---
    var p1 = document.createElement('p');
    p1.textContent = 'Transformers process all tokens in parallel (unlike RNNs). This is great for speed but means they have no inherent sense of order. Without positional information, the model would treat \u201CATGC\u201D and \u201CGCTA\u201D identically.';
    content.appendChild(p1);

    var p2 = document.createElement('p');
    p2.textContent = 'The original transformer uses sinusoidal positional encodings: fixed (not learned) signals added to each embedding. Each position gets a unique \u201Cfingerprint\u201D composed of sine and cosine waves at different frequencies. The different frequency bands operate like a binary encoding of position.';
    content.appendChild(p2);

    var p3 = document.createElement('p');
    p3.textContent = 'Some modern models use learned positional embeddings (BERT), rotary embeddings (RoPE in LLaMA), or ALiBi instead. The principle is the same: inject order information so the model knows where each token sits in the sequence.';
    content.appendChild(p3);

    // --- Biology analogy ---
    viz.renderAnalogyBox(
      'Position is everything in biology. The codon AUG at position 1 of an mRNA is a start signal; AUG at position 150 just codes for methionine. A cysteine at position 14 vs. position 88 determines which disulfide bridges form. Positional encoding captures this principle: same token, different position, different meaning.',
      content
    );

    // --- Sinusoidal canvas ---
    var canvasTitle = document.createElement('h3');
    canvasTitle.textContent = 'Sinusoidal Positional Encoding Waves';
    canvasTitle.style.cssText = 'margin-top: var(--space-xl); margin-bottom: var(--space-sm);';
    content.appendChild(canvasTitle);

    var canvasHint = document.createElement('p');
    canvasHint.style.cssText = 'font-size: 0.85rem; color: var(--text-muted); margin-bottom: var(--space-sm);';
    canvasHint.textContent = 'Each colored line represents one dimension of the positional encoding across all positions. Low-frequency dimensions (dim 0\u20131) change slowly; high-frequency dimensions (dim 6\u20137) oscillate rapidly.';
    content.appendChild(canvasHint);

    canvasContainer = document.createElement('div');
    content.appendChild(canvasContainer);

    // --- Before / After matrices ---
    var baTitle = document.createElement('h3');
    baTitle.textContent = 'Before & After: How Position Breaks Symmetry';
    baTitle.style.cssText = 'margin-top: var(--space-xl); margin-bottom: var(--space-sm);';
    content.appendChild(baTitle);

    var baHint = document.createElement('p');
    baHint.style.cssText = 'font-size: 0.85rem; color: var(--text-muted); margin-bottom: var(--space-sm);';
    baHint.textContent = 'Identical tokens (highlighted) have identical embedding rows before positional encoding. After adding PE, each position is unique \u2014 the model can now distinguish repeated tokens.';
    content.appendChild(baHint);

    beforeAfterContainer = document.createElement('div');
    content.appendChild(beforeAfterContainer);

    // --- Formula block ---
    formulaContainer = document.createElement('div');
    formulaContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(formulaContainer);

    // Render all
    renderCanvas();
    renderBeforeAfter();
    renderFormula();
  }

  window.TP.steps[3] = {
    init: function() {},
    render: function() {
      buildContent();
    }
  };
})();
