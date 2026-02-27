// step-1-tokenization.js — Step 1: Input & Tokenization
(function() {
  'use strict';
  var viz = window.TP.viz;
  var state = window.TP.state;

  var content = null;
  var inputField = null;
  var errorMsg = null;
  var tokenContainer = null;
  var formulaContainer = null;

  var DEFAULTS = { dna: 'ATGCGATC', protein: 'MAGLVKRD' };
  var VALID_CHARS = {
    dna: 'ATGC',
    protein: 'ACDEFGHIKLMNPQRSTVWY'
  };

  // --- Compute and store transformer results ---
  function recompute() {
    var seq = state.sequence;
    var mode = state.mode;
    if (!seq || seq.length < 4) return;
    state.result = TP.computeTransformer(seq, mode);
  }

  // --- Validate input sequence ---
  function validateInput(value) {
    var mode = state.mode;
    var valid = VALID_CHARS[mode];
    var upper = value.toUpperCase();
    var filtered = '';
    for (var i = 0; i < upper.length; i++) {
      if (valid.indexOf(upper[i]) >= 0) {
        filtered += upper[i];
      }
    }
    return filtered;
  }

  function onInputChange() {
    var raw = inputField.value;
    var filtered = validateInput(raw);

    // Show error if characters were filtered
    if (filtered !== raw.toUpperCase()) {
      inputField.classList.add('error');
      errorMsg.textContent = 'Only ' + (state.mode === 'dna' ? 'A, T, G, C' : 'standard amino acids') + ' allowed.';
      errorMsg.style.color = '#d9534f';
      setTimeout(function() { inputField.classList.remove('error'); }, 300);
    } else {
      errorMsg.textContent = '';
    }

    // Enforce length
    if (filtered.length > 10) {
      filtered = filtered.substring(0, 10);
      errorMsg.textContent = 'Maximum 10 tokens.';
      errorMsg.style.color = '#d9534f';
    }

    inputField.value = filtered;

    if (filtered.length < 4) {
      errorMsg.textContent = 'Need at least 4 tokens (currently ' + filtered.length + ').';
      errorMsg.style.color = 'var(--text-muted)';
      tokenContainer.innerHTML = '';
      return;
    }

    errorMsg.textContent = filtered.length + ' tokens (' + (state.mode === 'dna' ? '|V|=4' : '|V|=20') + ')';
    errorMsg.style.color = 'var(--text-muted)';

    state.sequence = filtered;
    recompute();
    renderTokens();
    renderFormula();
  }

  function renderTokens() {
    if (!tokenContainer) return;
    var seq = state.sequence;
    var mode = state.mode;
    if (!seq || seq.length < 4) {
      tokenContainer.innerHTML = '';
      return;
    }

    var tokens = seq.split('');
    var ids = TP.tokenize(seq, mode);
    var maskIdx = (mode === 'protein') ? Math.floor(seq.length / 2) : undefined;

    viz.renderTokenRow(tokens, tokenContainer, {
      ids: ids,
      maskIdx: maskIdx,
      mode: mode
    });
  }

  function renderFormula() {
    if (!formulaContainer) return;
    var seq = state.sequence;
    var mode = state.mode;
    if (!seq || seq.length < 4) return;

    var tokens = seq.split('');
    var ids = TP.tokenize(seq, mode);
    var vocabSize = mode === 'dna' ? 4 : 20;

    viz.renderFormulaBlock([
      'Input:    "' + seq + '"',
      'Tokens:   [' + tokens.map(function(t) { return '"' + t + '"'; }).join(', ') + ']',
      'Token IDs: [' + ids.join(', ') + ']',
      'Vocabulary V, |V| = ' + vocabSize + (mode === 'protein' ? ', [MASK] at position ' + Math.floor(seq.length / 2) : '')
    ], formulaContainer);
  }

  function buildContent() {
    content = document.getElementById('step-1-content');
    content.innerHTML = '';

    // --- Text ---
    var p1 = document.createElement('p');
    p1.textContent = 'Transformers operate on discrete tokens, not raw characters or continuous signals. Tokenization converts a raw input into a sequence of integer IDs from a fixed vocabulary. For DNA, our vocabulary is minimal: {A:0, T:1, G:2, C:3}. For proteins, we use the 20 standard amino acids.';
    content.appendChild(p1);

    var p2 = document.createElement('p');
    p2.textContent = 'Real genomic models sometimes use k-mer tokenization (e.g., DNABERT uses overlapping 6-mers), which is analogous to reading codons rather than individual bases. Real protein models include special tokens like [CLS] for classification, [SEP] for separation, [MASK] for masked language modeling, and [PAD] for padding shorter sequences.';
    content.appendChild(p2);

    // --- Biology analogy ---
    viz.renderAnalogyBox(
      'Tokenization is analogous to how the ribosome reads mRNA in discrete codon units. Just as the genetic code maps 64 codons to 20 amino acids (a many-to-one vocabulary mapping), tokenizers map raw sequences to a fixed set of IDs.',
      content
    );

    // --- Input field ---
    var inputCard = document.createElement('div');
    inputCard.className = 'tp-card';
    inputCard.style.textAlign = 'center';

    var inputLabel = document.createElement('div');
    inputLabel.style.cssText = 'font-family: var(--font-display); font-size: 0.82rem; font-weight: 600; color: var(--text-secondary); margin-bottom: var(--space-sm);';
    inputLabel.textContent = 'Enter a ' + (state.mode === 'dna' ? 'DNA' : 'protein') + ' sequence (4-10 characters):';
    inputLabel.id = 'step1-input-label';
    inputCard.appendChild(inputLabel);

    inputField = document.createElement('input');
    inputField.className = 'tp-input-field';
    inputField.type = 'text';
    inputField.maxLength = 10;
    inputField.value = state.sequence;
    inputField.placeholder = state.mode === 'dna' ? 'e.g. ATGCGATC' : 'e.g. MAGLVKRD';
    inputField.addEventListener('input', onInputChange);
    inputCard.appendChild(inputField);

    errorMsg = document.createElement('div');
    errorMsg.className = 'tp-input-help';
    errorMsg.textContent = state.sequence.length + ' tokens (' + (state.mode === 'dna' ? '|V|=4' : '|V|=20') + ')';
    inputCard.appendChild(errorMsg);

    content.appendChild(inputCard);

    // --- Token visualization ---
    var tokenTitle = document.createElement('h3');
    tokenTitle.textContent = 'Tokenized Sequence';
    tokenTitle.style.cssText = 'margin-top: var(--space-xl); margin-bottom: var(--space-sm);';
    content.appendChild(tokenTitle);

    tokenContainer = document.createElement('div');
    content.appendChild(tokenContainer);

    // --- Formula block ---
    formulaContainer = document.createElement('div');
    content.appendChild(formulaContainer);

    // --- Mode info ---
    var modeNote = document.createElement('p');
    modeNote.style.cssText = 'font-size: 0.85rem; color: var(--text-muted); margin-top: var(--space-lg);';
    if (state.mode === 'dna') {
      modeNote.innerHTML = '<strong>DNA mode (Decoder/GPT):</strong> Autoregressive — the model predicts the next nucleotide given all previous ones. A causal mask prevents each position from seeing future tokens.';
    } else {
      modeNote.innerHTML = '<strong>Protein mode (Encoder/BERT):</strong> Masked language modeling — one residue is hidden with [MASK] and the model predicts it from surrounding context. All positions can attend to all others.';
    }
    content.appendChild(modeNote);

    // Initial computation and render
    recompute();
    renderTokens();
    renderFormula();
  }

  window.TP.steps[1] = {
    init: function() {
      // Do initial computation so downstream steps have data
      recompute();
    },
    render: function() {
      buildContent();
    }
  };
})();
