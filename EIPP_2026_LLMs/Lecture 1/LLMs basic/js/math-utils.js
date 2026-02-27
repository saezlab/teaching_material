// math-utils.js — Linear algebra utilities, vocabulary maps, tokenization
// Pure functions, no DOM. All registered on window.TP.

(function() {
  'use strict';
  var TP = window.TP;

  // --- Toy model dimension constants ---
  TP.D_MODEL = 8;
  TP.N_HEADS = 4;
  TP.D_K = 2;       // d_model / n_heads
  TP.D_FF = 32;     // 4 * d_model

  // --- Vocabulary ---
  TP.DNA_VOCAB = { A: 0, T: 1, G: 2, C: 3 };
  TP.PROTEIN_VOCAB = {
    A:0, C:1, D:2, E:3, F:4, G:5, H:6, I:7, K:8, L:9,
    M:10, N:11, P:12, Q:13, R:14, S:15, T:16, V:17, W:18, Y:19
  };
  TP.DNA_VOCAB_REV = ['A','T','G','C'];
  TP.PROTEIN_VOCAB_REV = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'];

  // --- Tokenization ---
  TP.tokenize = function(seq, mode) {
    var vocab = mode === 'dna' ? TP.DNA_VOCAB : TP.PROTEIN_VOCAB;
    return seq.toUpperCase().split('').map(function(ch) {
      var id = vocab[ch];
      return id !== undefined ? id : -1;
    });
  };

  TP.detokenize = function(ids, mode) {
    var rev = mode === 'dna' ? TP.DNA_VOCAB_REV : TP.PROTEIN_VOCAB_REV;
    return ids.map(function(id) {
      return (id >= 0 && id < rev.length) ? rev[id] : '?';
    }).join('');
  };

  // --- Matrix / Vector Utilities ---

  // Create a 2D array of zeros: rows x cols
  TP.zeros = function(rows, cols) {
    var mat = [];
    for (var i = 0; i < rows; i++) {
      var row = [];
      for (var j = 0; j < cols; j++) {
        row.push(0);
      }
      mat.push(row);
    }
    return mat;
  };

  // Create a 1D array of zeros
  TP.zerosVec = function(n) {
    var arr = [];
    for (var i = 0; i < n; i++) arr.push(0);
    return arr;
  };

  // Matrix multiplication: A (m x p) * B (p x n) -> (m x n)
  TP.matMul = function(A, B) {
    var m = A.length;
    var p = A[0].length;
    var p2 = B.length;
    var n = B[0].length;
    if (p !== p2) {
      throw new Error('matMul dimension mismatch: A is ' + m + 'x' + p + ', B is ' + p2 + 'x' + n);
    }
    var C = TP.zeros(m, n);
    for (var i = 0; i < m; i++) {
      for (var j = 0; j < n; j++) {
        var sum = 0;
        for (var k = 0; k < p; k++) {
          sum += A[i][k] * B[k][j];
        }
        C[i][j] = sum;
      }
    }
    return C;
  };

  // Transpose a 2D array
  TP.transpose = function(A) {
    var rows = A.length;
    var cols = A[0].length;
    var T = TP.zeros(cols, rows);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < cols; j++) {
        T[j][i] = A[i][j];
      }
    }
    return T;
  };

  // Matrix element-wise addition: A + B (same dimensions)
  TP.matAdd = function(A, B) {
    var rows = A.length;
    var cols = A[0].length;
    var C = TP.zeros(rows, cols);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < cols; j++) {
        C[i][j] = A[i][j] + B[i][j];
      }
    }
    return C;
  };

  // Add bias vector to each row of a matrix: A (m x n) + b (length n)
  TP.addBias = function(A, b) {
    var rows = A.length;
    var cols = A[0].length;
    var C = TP.zeros(rows, cols);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < cols; j++) {
        C[i][j] = A[i][j] + b[j];
      }
    }
    return C;
  };

  // Dot product of two 1D arrays
  TP.dot = function(a, b) {
    var sum = 0;
    for (var i = 0; i < a.length; i++) {
      sum += a[i] * b[i];
    }
    return sum;
  };

  // Vector magnitude (L2 norm)
  TP.magnitude = function(a) {
    return Math.sqrt(TP.dot(a, a));
  };

  // Cosine similarity between two 1D arrays
  TP.cosineSimilarity = function(a, b) {
    var magA = TP.magnitude(a);
    var magB = TP.magnitude(b);
    if (magA === 0 || magB === 0) return 0;
    return TP.dot(a, b) / (magA * magB);
  };

  // Argmax of a 1D array — returns index of maximum value
  TP.argmax = function(arr) {
    var maxIdx = 0;
    var maxVal = arr[0];
    for (var i = 1; i < arr.length; i++) {
      if (arr[i] > maxVal) {
        maxVal = arr[i];
        maxIdx = i;
      }
    }
    return maxIdx;
  };

  // Scale a 1D or 2D array by a scalar
  TP.scale = function(A, s) {
    if (Array.isArray(A[0])) {
      // 2D
      return A.map(function(row) {
        return row.map(function(v) { return v * s; });
      });
    }
    // 1D
    return A.map(function(v) { return v * s; });
  };

  // --- Softmax (numerically stable) ---
  // Works on a 1D array. Subtracts max before exp to prevent overflow.
  TP.softmax = function(arr) {
    var max = arr[0];
    for (var i = 1; i < arr.length; i++) {
      if (arr[i] > max) max = arr[i];
    }
    var exps = [];
    var sum = 0;
    for (var i = 0; i < arr.length; i++) {
      var e = Math.exp(arr[i] - max);
      exps.push(e);
      sum += e;
    }
    return exps.map(function(e) { return e / sum; });
  };

  // --- ReLU ---
  // Works on 1D or 2D arrays
  TP.relu = function(A) {
    if (Array.isArray(A[0])) {
      // 2D
      return A.map(function(row) {
        return row.map(function(v) { return Math.max(0, v); });
      });
    }
    // 1D
    return A.map(function(v) { return Math.max(0, v); });
  };

  // --- Layer Normalization ---
  // Normalizes a 1D array (one token's representation across d_model dimensions)
  // gamma and beta are learnable scale/shift params (1D arrays of length d_model)
  // If not provided, defaults to gamma=1, beta=0 (standard normalization)
  TP.layerNorm = function(arr, gamma, beta, epsilon) {
    if (epsilon === undefined) epsilon = 1e-5;
    var n = arr.length;

    // Compute mean
    var sum = 0;
    for (var i = 0; i < n; i++) sum += arr[i];
    var mean = sum / n;

    // Compute variance
    var varSum = 0;
    for (var i = 0; i < n; i++) {
      var diff = arr[i] - mean;
      varSum += diff * diff;
    }
    var variance = varSum / n;

    // Normalize
    var result = [];
    var invStd = 1 / Math.sqrt(variance + epsilon);
    for (var i = 0; i < n; i++) {
      var normalized = (arr[i] - mean) * invStd;
      // Apply scale and shift if provided
      if (gamma && beta) {
        normalized = gamma[i] * normalized + beta[i];
      }
      result.push(normalized);
    }
    return result;
  };

  // Layer normalization for a 2D array (matrix): normalizes each row independently
  TP.layerNormMatrix = function(mat, gamma, beta, epsilon) {
    return mat.map(function(row) {
      return TP.layerNorm(row, gamma, beta, epsilon);
    });
  };

  // --- One-hot encoding ---
  TP.oneHot = function(id, vocabSize) {
    var vec = TP.zerosVec(vocabSize);
    if (id >= 0 && id < vocabSize) {
      vec[id] = 1;
    }
    return vec;
  };

  // --- Utility: get a row from a 2D array ---
  TP.getRow = function(mat, i) {
    return mat[i].slice();
  };

  // --- Utility: slice columns from a 2D array ---
  // Returns a new matrix with columns from colStart to colEnd (exclusive)
  TP.sliceCols = function(mat, colStart, colEnd) {
    return mat.map(function(row) {
      return row.slice(colStart, colEnd);
    });
  };

  // --- Utility: concatenate matrices horizontally ---
  // Takes an array of matrices with same number of rows, concatenates columns
  TP.concatCols = function(matrices) {
    var rows = matrices[0].length;
    var result = [];
    for (var i = 0; i < rows; i++) {
      var row = [];
      for (var m = 0; m < matrices.length; m++) {
        for (var j = 0; j < matrices[m][i].length; j++) {
          row.push(matrices[m][i][j]);
        }
      }
      result.push(row);
    }
    return result;
  };

  // --- Label helpers (shared across steps) ---
  TP.getDimLabels = function() {
    var labels = [];
    for (var i = 0; i < TP.D_MODEL; i++) {
      labels.push('d' + i);
    }
    return labels;
  };

  TP.getTokenLabels = function(tokens) {
    var labels = [];
    for (var i = 0; i < tokens.length; i++) {
      labels.push(i + ':' + tokens[i]);
    }
    return labels;
  };

})();
