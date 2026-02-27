// step-9-recap.js — Step 9: Recap & Real Models
(function() {
  'use strict';
  var viz = window.TP.viz;
  var state = window.TP.state;

  var FLOW_LABELS = [
    'Tokenization', 'Embedding', 'Positional Encoding', 'Self-Attention',
    'Multi-Head Attention', 'Residual + LN', 'Feed-Forward', 'Output'
  ];
  var FLOW_TO_STEP = [1, 2, 3, 4, 5, 6, 7, 8];

  var content = null;
  var flowContainer = null;
  var tableContainer = null;
  var formulaContainer = null;

  function renderFlowChart() {
    if (!flowContainer) return;
    flowContainer.innerHTML = '';

    var title = document.createElement('div');
    title.className = 'tp-viz-title';
    title.textContent = 'Complete Pipeline \u2014 click any box to revisit that step';
    flowContainer.appendChild(title);

    var diagramDiv = document.createElement('div');
    viz.renderFlowDiagram(FLOW_LABELS, 'all', diagramDiv, {
      onClick: function(idx) {
        if (TP.goToStep) {
          TP.goToStep(FLOW_TO_STEP[idx]);
        }
      }
    });
    flowContainer.appendChild(diagramDiv);
  }

  function renderTable() {
    if (!tableContainer) return;
    tableContainer.innerHTML = '';

    var title = document.createElement('h3');
    title.textContent = 'From Toy to Reality';
    title.style.cssText = 'margin-bottom: var(--space-sm);';
    tableContainer.appendChild(title);

    var hint = document.createElement('p');
    hint.style.cssText = 'font-size: 0.85rem; color: var(--text-muted); margin-bottom: var(--space-sm);';
    hint.textContent = 'Our toy model uses the exact same operations as real biological transformers \u2014 just at a much smaller scale.';
    tableContainer.appendChild(hint);

    // Table data
    var headers = ['', 'Toy model', 'ESM-2 (650M)', 'Nucleotide Transformer', 'DNABERT-2', 'AlphaFold Evoformer'];
    var rows = [
      ['d_model', '8', '1,280', '1,024', '768', '256'],
      ['Heads', '4', '20', '16', '12', '8'],
      ['Layers', '1', '33', '24', '12', '48'],
      ['d_ff', '32', '5,120', '4,096', '3,072', '1,024'],
      ['Params', '~500', '650M', '500M', '117M', '~90M'],
      ['Vocab', '4 / 20', '33', '4,096 (6-mers)', '4,096 (BPE)', '20 + struct']
    ];

    var table = document.createElement('table');
    table.className = 'tp-table';

    // Header row
    var thead = document.createElement('thead');
    var headerRow = document.createElement('tr');
    for (var h = 0; h < headers.length; h++) {
      var th = document.createElement('th');
      th.textContent = headers[h];
      if (h === 1) th.className = 'tp-table-highlight';
      headerRow.appendChild(th);
    }
    thead.appendChild(headerRow);
    table.appendChild(thead);

    // Body rows
    var tbody = document.createElement('tbody');
    for (var r = 0; r < rows.length; r++) {
      var tr = document.createElement('tr');
      for (var c = 0; c < rows[r].length; c++) {
        var td = document.createElement('td');
        td.textContent = rows[r][c];
        if (c === 0) td.style.fontWeight = '600';
        if (c === 1) td.className = 'tp-table-highlight';
        tr.appendChild(td);
      }
      tbody.appendChild(tr);
    }
    table.appendChild(tbody);

    var tableWrapper = document.createElement('div');
    tableWrapper.style.cssText = 'overflow-x: auto; border: 1px solid var(--card-border); border-radius: var(--border-radius);';
    tableWrapper.appendChild(table);
    tableContainer.appendChild(tableWrapper);
  }

  function renderFormula() {
    if (!formulaContainer) return;
    formulaContainer.innerHTML = '';

    viz.renderFormulaBlock([
      'Total parameters \u2248 L \u00D7 (4d\u00B2 + 2d\u00B7d_ff + 2d_ff + 6d)',
      'Where L = layers, d = d_model, d_ff = feed-forward dimension',
      'Toy: 1 \u00D7 (4\u00B78\u00B2 + 2\u00B78\u00B732 + 2\u00B732 + 6\u00B78) = 1 \u00D7 (256 + 512 + 64 + 48) \u2248 880',
      'ESM-2: 33 \u00D7 (4\u00D71280\u00B2 + ...) \u2248 650,000,000'
    ], formulaContainer);
  }

  function buildContent() {
    content = document.getElementById('step-9-content');
    content.innerHTML = '';

    // --- Summary text ---
    var p1 = document.createElement('p');
    p1.textContent = 'You\u2019ve now walked through every component of the transformer architecture. The key insight: the same sequence of operations \u2014 tokenization, embedding, positional encoding, self-attention, multi-head attention, residual connections, feed-forward networks, and output projection \u2014 powers both small toy models and state-of-the-art biological language models with hundreds of millions of parameters.';
    content.appendChild(p1);

    // --- Pipeline flow chart ---
    flowContainer = document.createElement('div');
    flowContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(flowContainer);

    // --- What changes at scale ---
    var h3Scale = document.createElement('h3');
    h3Scale.textContent = 'What Changes at Scale?';
    h3Scale.style.marginTop = 'var(--space-xl)';
    content.appendChild(h3Scale);

    var p2 = document.createElement('p');
    p2.textContent = 'More layers capture increasingly abstract representations. In protein models like ESM-2, lower layers capture local amino acid properties (charge, hydrophobicity); middle layers capture secondary structure (\u03B1-helices, \u03B2-sheets); upper layers capture tertiary contacts and functional sites. The same architecture, more compute, richer representations.';
    content.appendChild(p2);

    var p3 = document.createElement('p');
    p3.textContent = 'These models are pretrained on massive sequence databases (UniRef for proteins, GenBank for DNA) using self-supervised objectives \u2014 exactly the masked language modeling and next-token prediction you just explored. Fine-tuning adapts them to specific tasks: variant effect prediction, binding site detection, regulatory element classification, or protein structure prediction.';
    content.appendChild(p3);

    // --- Comparison table ---
    tableContainer = document.createElement('div');
    tableContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(tableContainer);

    // --- Biology analogy ---
    viz.renderAnalogyBox(
      'Our toy model is like a single E. coli protein: small, fast, and gets the point across. Real transformer models are like the human proteome \u2014 same fundamental building blocks (the same attention and FFN operations), but scaled up massively, enabling far more complex "functions." The leap from our 500-parameter toy to ESM-2\u2019s 650 million parameters mirrors the leap from a 300-amino-acid enzyme to the coordinated complexity of 20,000+ human proteins.',
      content
    );

    // --- Formula block ---
    formulaContainer = document.createElement('div');
    formulaContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(formulaContainer);

    // --- Closing note ---
    var closing = document.createElement('div');
    closing.style.cssText = 'margin-top: var(--space-xl); padding: var(--space-md); background: var(--accent-light); border-radius: var(--border-radius); text-align: center;';
    var closingText = document.createElement('p');
    closingText.style.cssText = 'margin: 0; font-family: var(--font-display); font-size: 0.95rem; font-weight: 600; color: var(--accent);';
    closingText.textContent = 'You now understand the transformer architecture from first principles. Every operation in GPT, BERT, ESM-2, and AlphaFold uses the same building blocks you explored here.';
    closing.appendChild(closingText);
    content.appendChild(closing);

    // Render all
    renderFlowChart();
    renderTable();
    renderFormula();
  }

  window.TP.steps[9] = {
    init: function() {},
    render: function() {
      buildContent();
    }
  };
})();
