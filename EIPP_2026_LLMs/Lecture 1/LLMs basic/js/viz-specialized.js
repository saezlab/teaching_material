// viz-specialized.js — Specialized visualization components
// Depends on viz-helpers.js.

(function() {
  'use strict';
  var viz = window.TP.viz;

  // =========================================================================
  // renderArrows — SVG attention arrows between token positions
  // =========================================================================
  // tokenRowEl: the .tp-token-row element returned by renderTokenRow
  // weights: array of attention weights [w0, w1, ...] (one per token)
  // sourceIdx: index of the source token (the row we're visualizing)
  // options: { threshold: 0.02 } — skip arrows below this weight
  viz.renderArrows = function(tokenRowEl, weights, sourceIdx, options) {
    options = options || {};
    var threshold = options.threshold || 0.02;

    // Remove any existing arrow SVG
    var parent = tokenRowEl.parentNode;
    var existing = parent.querySelector('.tp-arrows-svg');
    if (existing) existing.remove();

    // Get token block positions
    var tokens = tokenRowEl.querySelectorAll('.tp-token');
    if (tokens.length === 0) return;

    // Create SVG overlay
    var svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('class', 'tp-arrows-svg');

    // Position relative to parent
    var parentRect = parent.getBoundingClientRect();
    var rowRect = tokenRowEl.getBoundingClientRect();

    // SVG covers area above the token row
    var svgHeight = 80;
    svg.style.position = 'absolute';
    svg.style.left = (rowRect.left - parentRect.left) + 'px';
    svg.style.top = (rowRect.top - parentRect.top - svgHeight) + 'px';
    svg.style.width = rowRect.width + 'px';
    svg.style.height = svgHeight + 'px';
    svg.style.overflow = 'visible';
    svg.style.pointerEvents = 'none';
    svg.style.zIndex = '5';

    // Accent color from CSS
    var accent = getComputedStyle(document.documentElement).getPropertyValue('--accent').trim();

    for (var j = 0; j < weights.length; j++) {
      var w = weights[j];
      if (w < threshold) continue;

      var srcBlock = tokens[sourceIdx].querySelector('.tp-token-block');
      var tgtBlock = tokens[j].querySelector('.tp-token-block');
      if (!srcBlock || !tgtBlock) continue;

      var srcRect = srcBlock.getBoundingClientRect();
      var tgtRect = tgtBlock.getBoundingClientRect();

      // Coordinates relative to SVG
      var srcX = srcRect.left + srcRect.width / 2 - rowRect.left;
      var srcY = svgHeight; // bottom of svg = top of token row
      var tgtX = tgtRect.left + tgtRect.width / 2 - rowRect.left;
      var tgtY = svgHeight;

      // Draw curved arrow from source down to target
      var midY = svgHeight * (1 - Math.min(w * 1.5, 0.9));
      var opacity = 0.25 + w * 0.75;
      var strokeWidth = 1 + w * 5;

      var path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
      if (sourceIdx === j) {
        // Self-attention: small loop
        var d = 'M ' + srcX + ' ' + srcY +
                ' C ' + (srcX - 15) + ' ' + midY + ', ' +
                (srcX + 15) + ' ' + midY + ', ' +
                tgtX + ' ' + tgtY;
        path.setAttribute('d', d);
      } else {
        var d = 'M ' + srcX + ' ' + srcY +
                ' Q ' + ((srcX + tgtX) / 2) + ' ' + midY + ', ' +
                tgtX + ' ' + tgtY;
        path.setAttribute('d', d);
      }
      path.setAttribute('stroke', accent);
      path.setAttribute('stroke-width', strokeWidth);
      path.setAttribute('fill', 'none');
      path.setAttribute('opacity', opacity);
      path.setAttribute('stroke-linecap', 'round');

      // Arrowhead
      var arrow = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
      arrow.setAttribute('cx', tgtX);
      arrow.setAttribute('cy', tgtY - 2);
      arrow.setAttribute('r', 2 + w * 2);
      arrow.setAttribute('fill', accent);
      arrow.setAttribute('opacity', opacity);

      svg.appendChild(path);
      svg.appendChild(arrow);
    }

    // Make parent relative if not already
    if (getComputedStyle(parent).position === 'static') {
      parent.style.position = 'relative';
    }
    parent.appendChild(svg);
  };

  // =========================================================================
  // renderSinusoidalCanvas — Positional encoding wave visualization
  // =========================================================================
  // numPositions: number of positions (sequence length)
  // dModel: embedding dimension (8)
  // container: DOM element to append into
  // options: { width, height }
  viz.renderSinusoidalCanvas = function(numPositions, dModel, container, options) {
    options = options || {};
    container.innerHTML = '';

    var wrapper = document.createElement('div');
    wrapper.className = 'tp-canvas-container';

    var canvas = document.createElement('canvas');
    var width = options.width || 600;
    var height = options.height || 220;
    var dpr = window.devicePixelRatio || 1;
    canvas.width = width * dpr;
    canvas.height = height * dpr;
    canvas.style.width = width + 'px';
    canvas.style.height = height + 'px';

    var ctx = canvas.getContext('2d');
    ctx.scale(dpr, dpr);

    // Background
    ctx.fillStyle = '#fafaf8';
    ctx.fillRect(0, 0, width, height);

    // Margins
    var marginLeft = 50;
    var marginRight = 20;
    var marginTop = 20;
    var marginBottom = 40;
    var plotW = width - marginLeft - marginRight;
    var plotH = height - marginTop - marginBottom;

    // Axes
    ctx.strokeStyle = '#e8e6e1';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(marginLeft, marginTop);
    ctx.lineTo(marginLeft, marginTop + plotH);
    ctx.lineTo(marginLeft + plotW, marginTop + plotH);
    ctx.stroke();

    // Y-axis labels
    ctx.fillStyle = '#888';
    ctx.font = '11px JetBrains Mono, monospace';
    ctx.textAlign = 'right';
    ctx.fillText('1.0', marginLeft - 8, marginTop + 4);
    ctx.fillText('0.0', marginLeft - 8, marginTop + plotH / 2 + 4);
    ctx.fillText('-1.0', marginLeft - 8, marginTop + plotH + 4);

    // Zero line
    ctx.strokeStyle = '#e8e6e1';
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(marginLeft, marginTop + plotH / 2);
    ctx.lineTo(marginLeft + plotW, marginTop + plotH / 2);
    ctx.stroke();
    ctx.setLineDash([]);

    // X-axis labels (position indices)
    ctx.textAlign = 'center';
    ctx.fillStyle = '#888';
    for (var p = 0; p < numPositions; p++) {
      var x = marginLeft + (p / (numPositions - 1)) * plotW;
      ctx.fillText(p, x, marginTop + plotH + 18);
    }
    ctx.fillText('Position', marginLeft + plotW / 2, marginTop + plotH + 34);

    // Color palette for dimensions (8 distinct colors)
    var colors = [
      '#2d7d9a', '#e06c54', '#5ba355', '#c9873c',
      '#7b68ae', '#cb5a8c', '#3dacac', '#8e8e3e'
    ];

    // Draw PE curves
    for (var dim = 0; dim < dModel; dim++) {
      ctx.strokeStyle = colors[dim % colors.length];
      ctx.lineWidth = 2;
      ctx.beginPath();

      var nSamples = Math.max(numPositions * 10, 50);
      for (var s = 0; s <= nSamples; s++) {
        var pos = (s / nSamples) * (numPositions - 1);
        var i = Math.floor(dim / 2);
        var freq = 1.0 / Math.pow(10000, (2 * i) / dModel);
        var val;
        if (dim % 2 === 0) {
          val = Math.sin(pos * freq);
        } else {
          val = Math.cos(pos * freq);
        }

        var x = marginLeft + (pos / (numPositions - 1)) * plotW;
        var y = marginTop + plotH / 2 - val * (plotH / 2 - 4);

        if (s === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }

    wrapper.appendChild(canvas);

    // Legend
    var legend = document.createElement('div');
    legend.style.cssText = 'display:flex; flex-wrap:wrap; gap:12px; margin-top:8px; justify-content:center;';
    for (var dim = 0; dim < dModel; dim++) {
      var item = document.createElement('span');
      item.style.cssText = 'display:flex; align-items:center; gap:4px; font-family:var(--font-mono); font-size:0.7rem; color:#555;';
      var swatch = document.createElement('span');
      swatch.style.cssText = 'width:14px; height:3px; border-radius:2px; background:' + colors[dim % colors.length] + ';';
      item.appendChild(swatch);
      var funcType = dim % 2 === 0 ? 'sin' : 'cos';
      var freqIdx = Math.floor(dim / 2);
      item.appendChild(document.createTextNode('d' + dim + ' (' + funcType + ', i=' + freqIdx + ')'));
      legend.appendChild(item);
    }
    wrapper.appendChild(legend);

    container.appendChild(wrapper);
  };

  // =========================================================================
  // renderFlowDiagram — Vertical pipeline flow chart
  // =========================================================================
  // stepLabels: array of box labels
  // currentStep: index to highlight (-1 for none, 'all' to highlight all)
  // container: DOM element
  // options: { onClick: function(index) {}, highlightUpTo: number }
  viz.renderFlowDiagram = function(stepLabels, currentStep, container, options) {
    options = options || {};
    container.innerHTML = '';

    var flow = document.createElement('div');
    flow.className = 'tp-flow';

    for (var i = 0; i < stepLabels.length; i++) {
      if (i > 0) {
        // Arrow between boxes
        var arrow = document.createElement('div');
        arrow.className = 'tp-flow-arrow';
        flow.appendChild(arrow);
      }

      var box = document.createElement('div');
      box.className = 'tp-flow-box';
      box.textContent = stepLabels[i];

      // Highlighting logic
      if (currentStep === 'all') {
        box.classList.add('highlight');
      } else if (typeof options.highlightUpTo === 'number') {
        if (i <= options.highlightUpTo) {
          box.classList.add('highlight');
        }
      } else if (i === currentStep) {
        box.classList.add('highlight');
      }

      // Click handler
      if (options.onClick) {
        (function(idx) {
          box.addEventListener('click', function() {
            options.onClick(idx);
          });
        })(i);
      }

      flow.appendChild(box);
    }

    container.appendChild(flow);
  };

  // =========================================================================
  // renderNodeDiagram — FFN node visualization with play animation
  // =========================================================================
  // inputSize: number of input nodes (d_model = 8)
  // hiddenSubset: number of hidden nodes to show (representative ~12)
  // outputSize: number of output nodes (d_model = 8)
  // container: DOM element
  // options: {
  //   hiddenTotal: total hidden nodes (32),
  //   activations: array of hidden activations (to determine which are ReLU-zeroed),
  //   onPlay: callback when play animation is triggered
  // }
  viz.renderNodeDiagram = function(inputSize, hiddenSubset, outputSize, container, options) {
    options = options || {};
    container.innerHTML = '';

    var diagram = document.createElement('div');
    diagram.className = 'tp-node-diagram';

    var hiddenTotal = options.hiddenTotal || 32;
    var activations = options.activations || [];

    // Helper: create a column of nodes
    function makeCol(count, label, extraLabel) {
      var col = document.createElement('div');
      col.className = 'tp-node-col';
      var lbl = document.createElement('div');
      lbl.style.cssText = 'font-family:var(--font-display); font-size:0.72rem; font-weight:600; color:var(--text-muted); margin-bottom:4px;';
      lbl.textContent = label;
      col.appendChild(lbl);

      var nodes = [];
      for (var i = 0; i < count; i++) {
        var node = document.createElement('div');
        node.className = 'tp-node';
        col.appendChild(node);
        nodes.push(node);
      }

      if (extraLabel) {
        var extra = document.createElement('div');
        extra.style.cssText = 'font-family:var(--font-mono); font-size:0.7rem; color:var(--text-muted); margin-top:2px;';
        extra.textContent = extraLabel;
        col.appendChild(extra);
      }

      return { el: col, nodes: nodes };
    }

    // Connection lines container (SVG)
    var svgWrapper = document.createElement('div');
    svgWrapper.style.cssText = 'position:relative; display:flex; justify-content:center; align-items:center; gap:0;';

    var inputCol = makeCol(inputSize, 'Input', 'd=' + inputSize);
    var hiddenExtraLabel = hiddenSubset < hiddenTotal ? '... +' + (hiddenTotal - hiddenSubset) + ' more' : '';
    var hiddenCol = makeCol(hiddenSubset, 'Hidden (ReLU)', hiddenExtraLabel);
    var outputCol = makeCol(outputSize, 'Output', 'd=' + outputSize);

    // Arrow placeholders between columns
    var arrow1 = document.createElement('div');
    arrow1.style.cssText = 'font-size:1.4rem; color:var(--text-muted); padding:0 8px;';
    arrow1.textContent = '\u2192';
    var arrow2 = document.createElement('div');
    arrow2.style.cssText = 'font-size:1.4rem; color:var(--text-muted); padding:0 8px;';
    arrow2.textContent = '\u2192';

    diagram.appendChild(inputCol.el);
    diagram.appendChild(arrow1);
    diagram.appendChild(hiddenCol.el);
    diagram.appendChild(arrow2);
    diagram.appendChild(outputCol.el);

    container.appendChild(diagram);

    // Play button with sequential animation
    var playContainer = document.createElement('div');
    playContainer.style.cssText = 'text-align:center; margin-top:8px;';

    viz.createPlayButton(playContainer, function(done) {
      // Reset all nodes
      var allNodes = inputCol.nodes.concat(hiddenCol.nodes).concat(outputCol.nodes);
      allNodes.forEach(function(n) { n.classList.remove('lit', 'dark'); });

      // Phase 1: Light input nodes (staggered)
      var delay = 0;
      inputCol.nodes.forEach(function(n, i) {
        setTimeout(function() { n.classList.add('lit'); }, delay + i * 50);
      });
      delay += inputSize * 50 + 200;

      // Phase 2: Light hidden nodes (with ReLU zeroing)
      hiddenCol.nodes.forEach(function(n, i) {
        setTimeout(function() {
          // Determine if this node is active or ReLU-zeroed
          if (activations.length > 0) {
            // Sample from actual activations
            var actIdx = Math.floor(i * (activations.length / hiddenSubset));
            if (actIdx < activations.length && activations[actIdx] > 0) {
              n.classList.add('lit');
            } else {
              n.classList.add('dark');
            }
          } else {
            // Default: ~60% light, ~40% dark
            if (Math.random() > 0.4) {
              n.classList.add('lit');
            } else {
              n.classList.add('dark');
            }
          }
        }, delay + i * 40);
      });
      delay += hiddenSubset * 40 + 200;

      // Phase 3: Light output nodes
      outputCol.nodes.forEach(function(n, i) {
        setTimeout(function() { n.classList.add('lit'); }, delay + i * 50);
      });
      delay += outputSize * 50 + 100;

      setTimeout(done, delay);
    });

    container.appendChild(playContainer);
  };

  // =========================================================================
  // renderDistributionBars — Before/after inline bar charts
  // =========================================================================
  // before: array of values (pre-normalization)
  // after: array of values (post-normalization)
  // container: DOM element
  // options: { labels: ['d0','d1',...] }
  viz.renderDistributionBars = function(before, after, container, options) {
    options = options || {};
    container.innerHTML = '';

    var wrapper = document.createElement('div');
    wrapper.className = 'tp-dist-bars';

    function makeGroup(title, values, color) {
      var group = document.createElement('div');
      group.className = 'tp-dist-group';

      var h = document.createElement('h4');
      h.textContent = title;
      group.appendChild(h);

      // Find range for scaling
      var absMax = 0;
      for (var i = 0; i < values.length; i++) {
        if (Math.abs(values[i]) > absMax) absMax = Math.abs(values[i]);
      }
      if (absMax === 0) absMax = 1;

      var barContainer = document.createElement('div');
      barContainer.style.cssText = 'display:flex; gap:3px; align-items:center; height:100px; padding:4px 0;';

      for (var i = 0; i < values.length; i++) {
        var barWrapper = document.createElement('div');
        barWrapper.style.cssText = 'display:flex; flex-direction:column; align-items:center; flex:1; height:100%; justify-content:center; position:relative;';

        // Zero line reference
        var zeroLine = document.createElement('div');
        zeroLine.style.cssText = 'position:absolute; width:100%; height:1px; background:#ddd; top:50%;';
        barWrapper.appendChild(zeroLine);

        var bar = document.createElement('div');
        var normalizedH = Math.abs(values[i]) / absMax * 40;
        var isPositive = values[i] >= 0;

        bar.style.cssText = 'width:80%; max-width:20px; border-radius:2px; position:absolute; transition:height 0.4s; background:' + color + ';';
        bar.style.height = normalizedH + 'px';
        if (isPositive) {
          bar.style.bottom = '50%';
        } else {
          bar.style.top = '50%';
        }
        barWrapper.appendChild(bar);

        // Dim label
        var lbl = document.createElement('div');
        lbl.style.cssText = 'position:absolute; bottom:0; font-family:var(--font-mono); font-size:0.58rem; color:var(--text-muted);';
        lbl.textContent = options.labels ? options.labels[i] : i;
        barWrapper.appendChild(lbl);

        barContainer.appendChild(barWrapper);
      }

      group.appendChild(barContainer);

      // Stats
      var stats = document.createElement('div');
      stats.style.cssText = 'font-family:var(--font-mono); font-size:0.68rem; color:var(--text-muted); text-align:center; margin-top:4px;';
      var mean = 0, variance = 0;
      for (var i = 0; i < values.length; i++) mean += values[i];
      mean /= values.length;
      for (var i = 0; i < values.length; i++) variance += (values[i] - mean) * (values[i] - mean);
      variance /= values.length;
      stats.textContent = '\u03BC=' + mean.toFixed(2) + '  \u03C3\u00B2=' + variance.toFixed(2);
      group.appendChild(stats);

      return group;
    }

    var accent = getComputedStyle(document.documentElement).getPropertyValue('--accent').trim();
    wrapper.appendChild(makeGroup('Before (Residual Sum)', before, accent));
    wrapper.appendChild(makeGroup('After (Layer Norm)', after, '#5ba355'));

    container.appendChild(wrapper);
  };

})();
