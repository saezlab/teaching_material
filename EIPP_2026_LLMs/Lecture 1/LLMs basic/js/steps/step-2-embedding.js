// step-2-embedding.js — Step 2: Embedding Lookup
(function() {
  'use strict';
  var viz = window.TP.viz;
  var state = window.TP.state;

  var content = null;
  var fullMatrixContainer = null;
  var seqMatrixContainer = null;
  var formulaContainer = null;

  function getData() {
    if (!state.result) return null;
    return {
      seq: state.sequence,
      mode: state.mode,
      tokenIds: state.result.tokenIds,
      embeddings: state.result.embeddings
    };
  }

  function getFullEmbeddingMatrix() {
    var w = TP.WEIGHTS[state.mode];
    return w ? w.embedding : null;
  }

  function getVocabLabels() {
    if (state.mode === 'dna') {
      return ['A', 'T', 'G', 'C'];
    }
    // Protein: 20 amino acids in alphabetical order (matching TP.PROTEIN_VOCAB)
    return 'ACDEFGHIKLMNPQRSTVWY'.split('');
  }

  var getDimLabels = TP.getDimLabels;

  // Find which rows in the sequence have identical tokens
  function findIdenticalTokenRows(tokenIds) {
    var seen = {};
    var duplicateRows = [];
    for (var i = 0; i < tokenIds.length; i++) {
      var id = tokenIds[i];
      if (seen[id] !== undefined) {
        // Both the first occurrence and this one are duplicates
        if (duplicateRows.indexOf(seen[id]) < 0) {
          duplicateRows.push(seen[id]);
        }
        duplicateRows.push(i);
      } else {
        seen[id] = i;
      }
    }
    return duplicateRows;
  }

  function renderFullMatrix() {
    if (!fullMatrixContainer) return;
    var E = getFullEmbeddingMatrix();
    if (!E) return;

    var data = getData();
    // Highlight rows that correspond to tokens in the current sequence
    var highlightRows = [];
    if (data) {
      var unique = {};
      for (var i = 0; i < data.tokenIds.length; i++) {
        unique[data.tokenIds[i]] = true;
      }
      for (var id in unique) {
        highlightRows.push(parseInt(id));
      }
    }

    viz.renderMatrix(E, fullMatrixContainer, {
      title: 'Full Embedding Matrix E \u2208 \u211D' +
        (state.mode === 'dna' ? '\u2074' : '\u00B2\u2070') + '\u00D7\u2078',
      rowLabels: getVocabLabels(),
      colLabels: getDimLabels(),
      highlightRows: highlightRows,
      scale: 'blueWhiteRed'
    });
  }

  function renderSeqMatrix() {
    if (!seqMatrixContainer) return;
    var data = getData();
    if (!data || data.seq.length < 4) {
      seqMatrixContainer.innerHTML = '';
      return;
    }

    var tokens = data.seq.split('');
    var rowLabels = [];
    for (var i = 0; i < tokens.length; i++) {
      rowLabels.push(i + ':' + tokens[i]);
    }

    var duplicateRows = findIdenticalTokenRows(data.tokenIds);

    viz.renderMatrix(data.embeddings, seqMatrixContainer, {
      title: 'Sequence Embeddings (' + tokens.length + ' tokens \u00D7 ' + TP.D_MODEL + ' dims)',
      rowLabels: rowLabels,
      colLabels: getDimLabels(),
      highlightRows: duplicateRows,
      scale: 'blueWhiteRed'
    });

    // Annotation for identical tokens
    if (duplicateRows.length > 0) {
      var note = document.createElement('div');
      note.style.cssText = 'font-size: 0.8rem; color: var(--accent); margin-top: 6px; font-style: italic;';
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
          var positions = groups[id].join(', ');
          var token = tokens[groups[id][0]];
          parts.push('\u201C' + token + '\u201D at positions ' + positions);
        }
      }
      if (parts.length > 0) {
        note.textContent = 'Highlighted: identical tokens share identical embedding rows \u2014 ' + parts.join('; ') + '.';
        seqMatrixContainer.appendChild(note);
      }
    }
  }

  function renderFormula() {
    if (!formulaContainer) return;
    var data = getData();
    if (!data || data.seq.length < 4) return;

    var vocabSize = state.mode === 'dna' ? 4 : 20;

    viz.renderFormulaBlock([
      '<b>x</b><sub>embed</sub> = E[token_id] \u2208 \u211D<sup>d_model</sup>',
      'E \u2208 \u211D<sup>|V| \u00D7 d_model</sup> = \u211D<sup>' + vocabSize + ' \u00D7 ' + TP.D_MODEL + '</sup>',
      'Equivalent to: <b>x</b><sub>embed</sub> = onehot(token_id) \u00B7 E'
    ], formulaContainer);
  }

  function buildContent() {
    content = document.getElementById('step-2-content');
    content.innerHTML = '';

    // --- Text ---
    var p1 = document.createElement('p');
    p1.textContent = 'A token ID like "2" is meaningless to a neural network. The embedding layer converts each ID into a dense real-valued vector of dimension d_model. Mathematically, this is a lookup in a learned matrix E \u2208 \u211D^{|V| \u00D7 d_model}, or equivalently, multiplying the one-hot encoded token by E.';
    content.appendChild(p1);

    var p2 = document.createElement('p');
    p2.textContent = 'The key insight: embeddings place semantically or functionally similar tokens nearby in vector space. For proteins, we\u2019d expect hydrophobic residues to cluster and charged residues to cluster, similar to how BLOSUM62 captures substitution patterns from sequence alignments. Unlike BLOSUM (which is hand-crafted from MSA statistics), embeddings are learned end-to-end from the training objective and can capture subtler relationships.';
    content.appendChild(p2);

    // --- Biology analogy ---
    viz.renderAnalogyBox(
      'BLOSUM62 is a 20\u00D720 matrix encoding amino acid substitution probabilities, derived from blocks of aligned sequences. Embeddings are like a learned, higher-dimensional generalization: instead of 20\u00D720, we get 20\u00D7d_model, and the \u201Csimilarity\u201D is learned from billions of sequences rather than hand-curated alignments.',
      content
    );

    // --- Full embedding matrix ---
    var fullTitle = document.createElement('h3');
    fullTitle.textContent = 'Full Embedding Matrix';
    fullTitle.style.cssText = 'margin-top: var(--space-xl); margin-bottom: var(--space-sm);';
    content.appendChild(fullTitle);

    var fullHint = document.createElement('p');
    fullHint.style.cssText = 'font-size: 0.85rem; color: var(--text-muted); margin-bottom: var(--space-sm);';
    fullHint.textContent = 'Each row is one vocabulary token\u2019s embedding vector. Highlighted rows are tokens present in your current sequence. Hover any cell for exact values.';
    content.appendChild(fullHint);

    fullMatrixContainer = document.createElement('div');
    content.appendChild(fullMatrixContainer);

    // --- Sequence embeddings ---
    var seqTitle = document.createElement('h3');
    seqTitle.textContent = 'Sequence Embeddings';
    seqTitle.style.cssText = 'margin-top: var(--space-xl); margin-bottom: var(--space-sm);';
    content.appendChild(seqTitle);

    var seqHint = document.createElement('p');
    seqHint.style.cssText = 'font-size: 0.85rem; color: var(--text-muted); margin-bottom: var(--space-sm);';
    seqHint.textContent = 'The embedding vectors for each token in your input sequence. Note that identical tokens produce identical rows \u2014 the model cannot yet distinguish their positions.';
    content.appendChild(seqHint);

    seqMatrixContainer = document.createElement('div');
    content.appendChild(seqMatrixContainer);

    // --- Formula block ---
    formulaContainer = document.createElement('div');
    formulaContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(formulaContainer);

    // Render all
    renderFullMatrix();
    renderSeqMatrix();
    renderFormula();
  }

  window.TP.steps[2] = {
    init: function() {},
    render: function() {
      buildContent();
    }
  };
})();
