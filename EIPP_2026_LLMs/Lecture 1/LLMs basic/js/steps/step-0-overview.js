// step-0-overview.js — Step 0: Overview & Motivation
(function() {
  'use strict';
  var viz = window.TP.viz;
  var state = window.TP.state;

  var FLOW_LABELS = [
    'Tokenization', 'Embedding', 'Positional Encoding', 'Self-Attention',
    'Multi-Head Attention', 'Residual + LN', 'Feed-Forward', 'Output'
  ];

  // Map flow box index → step number (for click navigation)
  var FLOW_TO_STEP = [1, 2, 3, 4, 5, 6, 7, 8];

  var content = null;
  var flowContainer = null;

  function buildContent() {
    content = document.getElementById('step-0-content');
    content.innerHTML = '';

    // --- Paragraph 1: Long-range dependencies ---
    var p1 = document.createElement('p');
    p1.textContent = 'Biological sequences are full of long-range dependencies. A mutation at position 10 in a protein can disrupt folding at position 200. A regulatory element upstream of a gene can control expression thousands of bases away. In language, a pronoun at the end of a sentence refers back to a noun at the beginning. To model sequences like these, we need architectures that can capture relationships between distant positions efficiently.';
    content.appendChild(p1);

    // --- Paragraph 2: RNN → Attention history ---
    var p2 = document.createElement('p');
    p2.textContent = 'Early sequence models — Recurrent Neural Networks (RNNs) — processed sequences one element at a time, carrying forward a hidden state. This created an information bottleneck: by the time the model reaches position 200, it has largely "forgotten" position 10. LSTMs and GRUs improved on this with gating mechanisms, but still struggled with very long sequences. In 2017, the landmark paper "Attention Is All You Need" proposed the Transformer, which replaced recurrence entirely with a mechanism called self-attention — allowing every position to directly interact with every other position in a single step.';
    content.appendChild(p2);

    // --- Paragraph 3: Bioinformatics relevance ---
    var p3 = document.createElement('p');
    p3.textContent = 'This architecture now underpins the most powerful models in computational biology. Protein language models like ESM-2 (650M parameters) learn residue-level representations from millions of protein sequences. DNA language models like the Nucleotide Transformer and DNABERT-2 capture genomic patterns across thousands of bases. AlphaFold\'s Evoformer uses a specialized transformer to predict protein structure from sequence and evolutionary data. Understanding how transformers work is now essential for anyone working at the intersection of biology and machine learning.';
    content.appendChild(p3);

    // --- Biology analogy ---
    viz.renderAnalogyBox(
      'Think of an RNN like a ribosome: it reads the sequence one element at a time, carrying forward a hidden state. A transformer is more like a chaperone complex that can simultaneously consider all parts of the protein at once to determine the correct fold.',
      content
    );

    // --- Pipeline flow chart ---
    var flowTitle = document.createElement('h3');
    flowTitle.textContent = 'The Transformer Pipeline';
    flowTitle.style.cssText = 'margin-top: var(--space-xl); margin-bottom: var(--space-sm);';
    content.appendChild(flowTitle);

    var flowHint = document.createElement('p');
    flowHint.style.cssText = 'font-size: 0.85rem; color: var(--text-muted); margin-bottom: var(--space-sm);';
    flowHint.textContent = 'Click any box to jump to that step. Boxes highlight as you progress through the tutorial.';
    content.appendChild(flowHint);

    flowContainer = document.createElement('div');
    content.appendChild(flowContainer);

    renderFlow();
  }

  function renderFlow() {
    if (!flowContainer) return;

    // Highlight boxes up to the current step (step 0 = no boxes highlighted,
    // step 1 = box 0 highlighted, step 8 = all boxes highlighted)
    var highlightIdx = state.currentStep > 0 ? state.currentStep - 1 : -1;

    viz.renderFlowDiagram(FLOW_LABELS, highlightIdx, flowContainer, {
      highlightUpTo: highlightIdx >= 0 ? highlightIdx : undefined,
      onClick: function(idx) {
        TP.goToStep(FLOW_TO_STEP[idx]);
      }
    });
  }

  window.TP.steps[0] = {
    init: function() {
      buildContent();
    },
    render: function() {
      // Re-render flow chart to update highlights based on progress
      if (!flowContainer) {
        buildContent();
      } else {
        renderFlow();
      }
    }
  };
})();
