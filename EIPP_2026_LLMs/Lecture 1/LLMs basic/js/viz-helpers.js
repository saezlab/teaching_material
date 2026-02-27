// viz-helpers.js — Core reusable visualization renderers
// Registered on TP.viz. CSS lives in transformer_playground.html <style>.

(function() {
  'use strict';
  var viz = window.TP.viz;

  // =========================================================================
  // Shared tooltip element (one global tooltip, repositioned on hover)
  // =========================================================================
  var tooltip = null;
  function getTooltip() {
    if (!tooltip) {
      tooltip = document.createElement('div');
      tooltip.className = 'tp-tooltip';
      document.body.appendChild(tooltip);
    }
    return tooltip;
  }

  function showTooltip(text, e) {
    var tt = getTooltip();
    tt.textContent = text;
    tt.classList.add('visible');
    tt.style.left = (e.clientX + 12) + 'px';
    tt.style.top = (e.clientY - 28) + 'px';
  }

  function hideTooltip() {
    var tt = getTooltip();
    tt.classList.remove('visible');
  }

  // =========================================================================
  // Color mapping functions
  // =========================================================================

  // Blue (negative) → White (zero) → Red (positive)
  function valueToBlueWhiteRed(v, min, max) {
    // Normalize to [-1, 1]
    var range = Math.max(Math.abs(min), Math.abs(max));
    if (range === 0) return 'rgb(255,255,255)';
    var t = v / range; // -1 to 1

    if (t < 0) {
      // Blue to white: t goes from -1 (full blue) to 0 (white)
      var s = -t; // 0 to 1
      var r = Math.round(255 * (1 - s * 0.7));
      var g = Math.round(255 * (1 - s * 0.7));
      var b = 255;
      return 'rgb(' + r + ',' + g + ',' + b + ')';
    } else {
      // White to red: t goes from 0 (white) to 1 (full red)
      var s = t;
      var r = 255;
      var g = Math.round(255 * (1 - s * 0.7));
      var b = Math.round(255 * (1 - s * 0.7));
      return 'rgb(' + r + ',' + g + ',' + b + ')';
    }
  }

  // Cream → Green (5-stop scale for attention weights 0→1)
  // Stops: 0=#faf6e8 (cream), 0.25=#d4e8b0, 0.5=#8cc379, 0.75=#4a9f55, 1=#1a6d32
  function valueToGreenScale(v) {
    v = Math.max(0, Math.min(1, v));
    var stops = [
      [250, 246, 232],  // cream
      [212, 232, 176],  // light green
      [140, 195, 121],  // medium green
      [74,  159, 85],   // dark green
      [26,  109, 50]    // very dark green
    ];
    var t = v * 4;
    var i = Math.min(3, Math.floor(t));
    var f = t - i;
    var c0 = stops[i];
    var c1 = stops[i + 1];
    var r = Math.round(c0[0] + (c1[0] - c0[0]) * f);
    var g = Math.round(c0[1] + (c1[1] - c0[1]) * f);
    var b = Math.round(c0[2] + (c1[2] - c0[2]) * f);
    return 'rgb(' + r + ',' + g + ',' + b + ')';
  }

  // Text color for readability on colored backgrounds
  function textColorForBg(v, isAttention) {
    if (isAttention) {
      return v > 0.55 ? '#fff' : '#333';
    }
    return (Math.abs(v) > 0.6) ? '#fff' : '#333';
  }

  // Expose color functions for use in specialized viz
  viz.valueToBlueWhiteRed = valueToBlueWhiteRed;
  viz.valueToGreenScale = valueToGreenScale;

  // =========================================================================
  // renderMatrix — Color-coded matrix / heatmap
  // =========================================================================
  // options: {
  //   scale: 'blueWhiteRed' (default) | 'greenScale',
  //   rowLabels: ['A','T',...],
  //   colLabels: ['d0','d1',...],
  //   onClick: function(row, col, value) {},
  //   title: 'Matrix title',
  //   highlightRows: [2, 4]  // row indices to add highlight border
  // }
  viz.renderMatrix = function(data, container, options) {
    options = options || {};
    container.innerHTML = '';

    var rows = data.length;
    var cols = data[0].length;
    var hasRowLabels = options.rowLabels && options.rowLabels.length > 0;
    var hasColLabels = options.colLabels && options.colLabels.length > 0;
    var isAttention = options.scale === 'greenScale';

    // Find min/max for blue-white-red scale
    var min = 0, max = 0;
    if (!isAttention) {
      for (var i = 0; i < rows; i++) {
        for (var j = 0; j < cols; j++) {
          if (data[i][j] < min) min = data[i][j];
          if (data[i][j] > max) max = data[i][j];
        }
      }
    }

    // Title
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

    var gridCols = cols + (hasRowLabels ? 1 : 0);
    var gridRows = rows + (hasColLabels ? 1 : 0);
    grid.style.gridTemplateColumns = 'repeat(' + gridCols + ', auto)';

    // Column labels row
    if (hasColLabels) {
      if (hasRowLabels) {
        var corner = document.createElement('div');
        corner.className = 'tp-matrix-label';
        grid.appendChild(corner);
      }
      for (var j = 0; j < cols; j++) {
        var cl = document.createElement('div');
        cl.className = 'tp-matrix-label';
        cl.textContent = options.colLabels[j];
        grid.appendChild(cl);
      }
    }

    // Data rows
    for (var i = 0; i < rows; i++) {
      if (hasRowLabels) {
        var rl = document.createElement('div');
        rl.className = 'tp-matrix-label';
        rl.textContent = options.rowLabels[i];
        if (options.highlightRows && options.highlightRows.indexOf(i) >= 0) {
          rl.style.fontWeight = '700';
          rl.style.color = 'var(--accent)';
        }
        grid.appendChild(rl);
      }

      for (var j = 0; j < cols; j++) {
        var cell = document.createElement('div');
        cell.className = 'tp-matrix-cell';
        var v = data[i][j];

        if (isAttention) {
          cell.style.backgroundColor = valueToGreenScale(v);
          cell.style.color = textColorForBg(v, true);
        } else {
          cell.style.backgroundColor = valueToBlueWhiteRed(v, min, max);
          cell.style.color = textColorForBg(v, false);
        }

        cell.textContent = v.toFixed(2);

        if (options.highlightRows && options.highlightRows.indexOf(i) >= 0) {
          cell.style.outline = '2px solid var(--accent)';
          cell.style.outlineOffset = '-2px';
        }

        // Tooltip
        (function(row, col, val) {
          cell.addEventListener('mouseenter', function(e) {
            var label = '';
            if (options.rowLabels) label += options.rowLabels[row];
            if (options.colLabels) label += (label ? ',' : '') + options.colLabels[col];
            showTooltip((label ? '[' + label + '] ' : '') + val.toFixed(4), e);
          });
          cell.addEventListener('mousemove', function(e) {
            var tt = getTooltip();
            tt.style.left = (e.clientX + 12) + 'px';
            tt.style.top = (e.clientY - 28) + 'px';
          });
          cell.addEventListener('mouseleave', hideTooltip);

          if (options.onClick) {
            cell.style.cursor = 'pointer';
            cell.addEventListener('click', function() {
              options.onClick(row, col, val);
            });
          }
        })(i, j, v);

        grid.appendChild(cell);
      }
    }

    wrapper.appendChild(grid);
    container.appendChild(wrapper);
  };

  // =========================================================================
  // renderTokenRow — Colored token blocks
  // =========================================================================
  // options: {
  //   ids: [0,1,2,3],        // token IDs to display below
  //   maskIdx: 4,             // index to show [MASK] overlay
  //   mode: 'dna' | 'protein'
  // }
  viz.renderTokenRow = function(tokens, container, options) {
    options = options || {};
    container.innerHTML = '';

    var row = document.createElement('div');
    row.className = 'tp-token-row';

    // Generate slightly varied accent colors for visual interest
    var mode = options.mode || 'dna';

    for (var i = 0; i < tokens.length; i++) {
      var token = document.createElement('div');
      token.className = 'tp-token';

      var block = document.createElement('div');
      block.className = 'tp-token-block';

      var isMasked = (options.maskIdx !== undefined && i === options.maskIdx);

      if (isMasked) {
        block.classList.add('masked');
        // The ::after pseudo-element shows [M]
        block.textContent = '';
      } else {
        block.textContent = tokens[i];
      }

      // Vary the shade slightly per position for visual interest
      var hueShift = i * 8;
      if (mode === 'dna') {
        block.style.backgroundColor = 'hsl(' + (192 + hueShift % 20) + ', 51%, ' + (38 + (i % 3) * 5) + '%)';
      } else {
        block.style.backgroundColor = 'hsl(' + (25 + hueShift % 15) + ', 40%, ' + (38 + (i % 3) * 5) + '%)';
      }

      token.appendChild(block);

      // Token ID below
      if (options.ids) {
        var idEl = document.createElement('div');
        idEl.className = 'tp-token-id';
        idEl.textContent = isMasked ? '?' : options.ids[i];
        token.appendChild(idEl);
      }

      row.appendChild(token);
    }

    container.appendChild(row);
    return row;
  };

  // =========================================================================
  // renderBarChart — Flexbox bar chart with animated height
  // =========================================================================
  // options: {
  //   highlightIdx: 2,        // index to highlight
  //   maxHeight: 150,          // max bar height in px (default 150)
  //   animate: true            // animate height (default true)
  // }
  viz.renderBarChart = function(labels, values, container, options) {
    options = options || {};
    container.innerHTML = '';

    var chart = document.createElement('div');
    chart.className = 'tp-bar-chart';

    var maxVal = 0;
    for (var i = 0; i < values.length; i++) {
      if (values[i] > maxVal) maxVal = values[i];
    }
    if (maxVal === 0) maxVal = 1;

    var maxHeight = options.maxHeight || 150;

    for (var i = 0; i < values.length; i++) {
      var wrapper = document.createElement('div');
      wrapper.className = 'tp-bar-wrapper';

      // Value label above bar
      var valLabel = document.createElement('div');
      valLabel.style.cssText = 'font-family: var(--font-mono); font-size: 0.62rem; color: var(--text-muted); margin-bottom: 2px;';
      valLabel.textContent = values[i].toFixed(2);
      wrapper.appendChild(valLabel);

      var bar = document.createElement('div');
      bar.className = 'tp-bar';
      if (options.highlightIdx !== undefined && i === options.highlightIdx) {
        bar.classList.add('highlight');
      }

      var targetHeight = Math.max(2, (values[i] / maxVal) * maxHeight);

      if (options.animate !== false) {
        bar.style.height = '2px';
        (function(b, h) {
          setTimeout(function() { b.style.height = h + 'px'; }, 50 + i * 40);
        })(bar, targetHeight);
      } else {
        bar.style.height = targetHeight + 'px';
      }

      // Tooltip on hover
      (function(idx, val) {
        bar.addEventListener('mouseenter', function(e) {
          showTooltip(labels[idx] + ': ' + val.toFixed(4), e);
        });
        bar.addEventListener('mousemove', function(e) {
          var tt = getTooltip();
          tt.style.left = (e.clientX + 12) + 'px';
          tt.style.top = (e.clientY - 28) + 'px';
        });
        bar.addEventListener('mouseleave', hideTooltip);
      })(i, values[i]);

      wrapper.appendChild(bar);

      var label = document.createElement('div');
      label.className = 'tp-bar-label';
      label.textContent = labels[i];
      wrapper.appendChild(label);

      chart.appendChild(wrapper);
    }

    container.appendChild(chart);
  };

  // =========================================================================
  // renderAnalogyBox — Biology analogy callout
  // =========================================================================
  viz.renderAnalogyBox = function(text, container) {
    var box = document.createElement('div');
    box.className = 'tp-analogy';
    box.textContent = text;
    container.appendChild(box);
    return box;
  };

  // =========================================================================
  // renderFormulaBlock — JetBrains Mono formula display
  // =========================================================================
  // lines: array of HTML strings (supports <sub>, <sup>, <b>, <i>)
  viz.renderFormulaBlock = function(lines, container) {
    var block = document.createElement('div');
    block.className = 'tp-formula';

    for (var i = 0; i < lines.length; i++) {
      var line = document.createElement('div');
      line.innerHTML = lines[i];
      block.appendChild(line);
    }

    container.appendChild(block);
    return block;
  };

  // =========================================================================
  // createPlayButton — ▶ button that disables during playback
  // =========================================================================
  // onPlay: function(done) — called when clicked. Call done() when animation finishes.
  viz.createPlayButton = function(container, onPlay) {
    var btn = document.createElement('button');
    btn.className = 'tp-play-btn';
    btn.innerHTML = '&#9654; Play';

    btn.addEventListener('click', function() {
      btn.disabled = true;
      btn.innerHTML = '&#9654; Playing...';
      onPlay(function() {
        btn.disabled = false;
        btn.innerHTML = '&#9654; Replay';
      });
    });

    container.appendChild(btn);
    return btn;
  };

})();
