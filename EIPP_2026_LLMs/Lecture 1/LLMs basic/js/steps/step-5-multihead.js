// step-5-multihead.js — Step 5: Multi-Head Attention
(function() {
  'use strict';
  var viz = window.TP.viz;
  var state = window.TP.state;

  var content = null;
  var tabContainer = null;
  var heatmapContainer = null;
  var tokenRowContainer = null;
  var concatContainer = null;
  var formulaContainer = null;
  var tokenRowEl = null;
  var activeHead = 0;

  var HEAD_LABELS = {
    dna: [
      'Complementary base pairing',
      'Local/sequential context',
      'Purine-pyrimidine patterns',
      'Diffuse/distributed'
    ],
    protein: [
      'Hydrophobic clustering',
      'Charge interactions',
      'Local backbone context',
      'Diffuse/distributed'
    ]
  };

  function getData() {
    if (!state.result || !state.result.qkvPerHead) return null;
    return {
      seq: state.sequence,
      mode: state.mode,
      tokens: state.sequence.split(''),
      tokenIds: state.result.tokenIds,
      qkvPerHead: state.result.qkvPerHead
    };
  }

  var getTokenLabels = TP.getTokenLabels;

  function renderTabs() {
    if (!tabContainer) return;
    tabContainer.innerHTML = '';

    var tabs = document.createElement('div');
    tabs.className = 'tp-tabs';

    for (var h = 0; h < TP.N_HEADS; h++) {
      var tab = document.createElement('button');
      tab.className = 'tp-tab';
      if (h === activeHead) tab.classList.add('active');
      tab.textContent = 'Head ' + (h + 1);
      (function(idx) {
        tab.addEventListener('click', function() {
          activeHead = idx;
          renderTabs();
          renderHeatmap();
          renderTokenRow();
        });
      })(h);
      tabs.appendChild(tab);
    }

    tabContainer.appendChild(tabs);

    // Head description label
    var labels = HEAD_LABELS[state.mode] || HEAD_LABELS.dna;
    var desc = document.createElement('div');
    desc.style.cssText = 'font-size: 0.82rem; color: var(--accent); font-style: italic; margin-top: var(--space-xs); text-align: center;';
    desc.textContent = 'Head ' + (activeHead + 1) + ': ' + labels[activeHead];
    tabContainer.appendChild(desc);
  }

  function renderHeatmap() {
    if (!heatmapContainer) return;
    var data = getData();
    if (!data) { heatmapContainer.innerHTML = ''; return; }

    var attn = data.qkvPerHead[activeHead].attention;
    var tokenLabels = getTokenLabels(data.tokens);

    viz.renderMatrix(attn, heatmapContainer, {
      title: 'Attention Weights \u2014 Head ' + (activeHead + 1) + ' (click a cell for arrows)',
      rowLabels: tokenLabels,
      colLabels: tokenLabels,
      scale: 'greenScale',
      onClick: function(row, col, val) {
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

    var arrowWrapper = document.createElement('div');
    arrowWrapper.style.cssText = 'position: relative; padding-top: 90px;';

    var maskIdx = (data.mode === 'protein') ? Math.floor(data.seq.length / 2) : undefined;
    tokenRowEl = viz.renderTokenRow(data.tokens, arrowWrapper, {
      ids: data.tokenIds,
      maskIdx: maskIdx,
      mode: data.mode
    });

    tokenRowContainer.appendChild(arrowWrapper);
  }

  function renderConcatDiagram() {
    if (!concatContainer) return;
    concatContainer.innerHTML = '';

    var title = document.createElement('div');
    title.className = 'tp-viz-title';
    title.textContent = 'Concatenation + Output Projection';
    concatContainer.appendChild(title);

    // Visual: 4 colored blocks → concat bar → W_O → output bar
    var diagram = document.createElement('div');
    diagram.style.cssText = 'display: flex; align-items: center; gap: 8px; flex-wrap: wrap; justify-content: center; padding: var(--space-sm) 0;';

    var headColors = ['#4a9f55', '#2d7d9a', '#c27a3a', '#8b6bb5'];
    var headLabels = ['H1', 'H2', 'H3', 'H4'];

    // 4 head blocks
    for (var h = 0; h < 4; h++) {
      var block = document.createElement('div');
      block.style.cssText = 'width: 48px; height: 36px; background: ' + headColors[h] + '; color: #fff; border-radius: 4px; display: flex; align-items: center; justify-content: center; font-family: var(--font-mono); font-size: 0.7rem; font-weight: 600;';
      block.textContent = headLabels[h];
      diagram.appendChild(block);
    }

    // Arrow: →
    var arrow1 = document.createElement('div');
    arrow1.style.cssText = 'font-size: 1.2rem; color: var(--text-muted);';
    arrow1.textContent = '\u2192';
    diagram.appendChild(arrow1);

    // Concat bar (striped with head colors)
    var concatBar = document.createElement('div');
    concatBar.style.cssText = 'width: 200px; height: 36px; border-radius: 4px; display: flex; overflow: hidden;';
    for (var h = 0; h < 4; h++) {
      var seg = document.createElement('div');
      seg.style.cssText = 'flex: 1; background: ' + headColors[h] + '; opacity: 0.7;';
      concatBar.appendChild(seg);
    }
    diagram.appendChild(concatBar);

    // Arrow: → W_O →
    var arrow2 = document.createElement('div');
    arrow2.style.cssText = 'font-size: 0.75rem; color: var(--text-muted); font-family: var(--font-mono);';
    arrow2.innerHTML = '\u00D7 W<sub>O</sub> \u2192';
    diagram.appendChild(arrow2);

    // Output bar
    var outBar = document.createElement('div');
    outBar.style.cssText = 'width: 80px; height: 36px; background: var(--accent); opacity: 0.8; border-radius: 4px; display: flex; align-items: center; justify-content: center; color: #fff; font-family: var(--font-mono); font-size: 0.7rem; font-weight: 600;';
    outBar.textContent = 'd_model';
    diagram.appendChild(outBar);

    concatContainer.appendChild(diagram);

    // Dimension annotation
    var dims = document.createElement('div');
    dims.style.cssText = 'font-size: 0.75rem; color: var(--text-muted); text-align: center; margin-top: 4px; font-family: var(--font-mono);';
    dims.textContent = '4 \u00D7 d_k=' + TP.D_K + ' = ' + (4 * TP.D_K) + ' \u2192 W_O \u2192 d_model=' + TP.D_MODEL;
    concatContainer.appendChild(dims);
  }

  function renderFormula() {
    if (!formulaContainer) return;
    formulaContainer.innerHTML = '';

    viz.renderFormulaBlock([
      'head<sub>i</sub> = Attention(XW<sub>Q</sub><sup>i</sup>, XW<sub>K</sub><sup>i</sup>, XW<sub>V</sub><sup>i</sup>)  for i = 1..' + TP.N_HEADS,
      'MultiHead(X) = Concat(head<sub>1</sub>, ..., head<sub>' + TP.N_HEADS + '</sub>) \u00B7 W<sub>O</sub>',
      'W<sub>O</sub> \u2208 \u211D<sup>' + (TP.N_HEADS * TP.D_K) + ' \u00D7 ' + TP.D_MODEL + '</sup>  (h\u00B7d_k \u00D7 d_model)',
      'Total parameters per multi-head attention: 4 \u00B7 d_model\u00B2 = ' + (4 * TP.D_MODEL * TP.D_MODEL)
    ], formulaContainer);
  }

  function buildContent() {
    content = document.getElementById('step-5-content');
    content.innerHTML = '';
    activeHead = 0;

    // --- Text ---
    var p1 = document.createElement('p');
    p1.textContent = 'A single attention head can only focus on one type of relationship per layer. Multi-head attention runs ' + TP.N_HEADS + ' parallel attention operations, each with its own learned W_Q, W_K, W_V projections. This allows the model to simultaneously attend to information from different \u201Csubspaces.\u201D';
    content.appendChild(p1);

    var p2 = document.createElement('p');
    p2.textContent = 'Each head operates on a reduced dimension d_k = d_model / h = ' + TP.D_K + ', so the total computation is the same as single-head attention with full d_model. The outputs are concatenated and projected back to d_model via W_O. In practice, different heads empirically specialize: in protein models like ESM-2, researchers have shown that attention heads correlate with contact maps, binding sites, and secondary structure.';
    content.appendChild(p2);

    // --- Biology analogy ---
    viz.renderAnalogyBox(
      'Proteins are stabilized by multiple types of non-covalent interactions simultaneously: hydrogen bonds, van der Waals forces, hydrophobic packing, electrostatic attractions, and \u03C0-stacking. No single \u201Cinteraction function\u201D captures all of these. Multi-head attention is analogous: each head specializes in detecting a different type of relationship.',
      content
    );

    // --- Tabbed interface ---
    var tabCard = document.createElement('div');
    tabCard.className = 'tp-card';
    tabCard.style.marginTop = 'var(--space-lg)';

    tabContainer = document.createElement('div');
    tabCard.appendChild(tabContainer);

    heatmapContainer = document.createElement('div');
    tabCard.appendChild(heatmapContainer);

    content.appendChild(tabCard);

    // --- Token row for arrows ---
    tokenRowContainer = document.createElement('div');
    content.appendChild(tokenRowContainer);

    // --- Concatenation diagram ---
    concatContainer = document.createElement('div');
    concatContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(concatContainer);

    // --- Formula block ---
    formulaContainer = document.createElement('div');
    formulaContainer.style.marginTop = 'var(--space-lg)';
    content.appendChild(formulaContainer);

    // Render all
    renderTabs();
    renderHeatmap();
    renderTokenRow();
    renderConcatDiagram();
    renderFormula();
  }

  window.TP.steps[5] = {
    init: function() {},
    render: function() {
      buildContent();
    }
  };
})();
