// transformer.js — Full transformer forward pass
// Depends on math-utils.js and weights.js.
// Returns ALL intermediate representations for visualization.

(function() {
  'use strict';
  var TP = window.TP;

  // --- Sinusoidal Positional Encoding ---
  // PE(pos, 2i) = sin(pos / 10000^(2i/d_model))
  // PE(pos, 2i+1) = cos(pos / 10000^(2i/d_model))
  function computePositionalEncoding(seqLen, dModel) {
    var PE = TP.zeros(seqLen, dModel);
    for (var pos = 0; pos < seqLen; pos++) {
      for (var i = 0; i < dModel; i += 2) {
        var angle = pos / Math.pow(10000, i / dModel);
        PE[pos][i] = Math.sin(angle);
        if (i + 1 < dModel) {
          PE[pos][i + 1] = Math.cos(angle);
        }
      }
    }
    return PE;
  }

  // --- Single-Head Attention ---
  // Returns { Q, K, V, scores, scaledScores, mask, attention, headOutput }
  function computeHeadAttention(X, headWeights, seqLen, dK, causal) {
    // Q = X * W_Q, K = X * W_K, V = X * W_V
    var Q = TP.matMul(X, headWeights.W_Q);
    var K = TP.matMul(X, headWeights.W_K);
    var V = TP.matMul(X, headWeights.W_V);

    // scores = Q * K^T
    var scores = TP.matMul(Q, TP.transpose(K));

    // scale by 1/sqrt(d_k)
    var scale = 1.0 / Math.sqrt(dK);
    var scaledScores = TP.scale(scores, scale);

    // Apply causal mask (DNA mode): set future positions to -Infinity
    var mask = TP.zeros(seqLen, seqLen);
    var maskedScores = [];
    for (var i = 0; i < seqLen; i++) {
      var row = [];
      for (var j = 0; j < seqLen; j++) {
        if (causal && j > i) {
          mask[i][j] = 1; // 1 = masked
          row.push(-Infinity);
        } else {
          row.push(scaledScores[i][j]);
        }
      }
      maskedScores.push(row);
    }

    // Softmax per row → attention weights
    var attention = [];
    for (var i = 0; i < seqLen; i++) {
      attention.push(TP.softmax(maskedScores[i]));
    }

    // output = attention * V
    var headOutput = TP.matMul(attention, V);

    return {
      Q: Q,
      K: K,
      V: V,
      scores: scores,
      scaledScores: scaledScores,
      mask: mask,
      attention: attention,
      headOutput: headOutput
    };
  }

  // --- Full Transformer Forward Pass ---
  TP.computeTransformer = function(sequence, mode) {
    var weights = TP.WEIGHTS[mode];
    var dModel = TP.D_MODEL;
    var nHeads = TP.N_HEADS;
    var dK = TP.D_K;
    var causal = (mode === 'dna');

    // Step 1: Tokenization
    var tokenIds = TP.tokenize(sequence, mode);
    var seqLen = tokenIds.length;

    // Step 2: Embedding lookup
    var embeddings = [];
    for (var i = 0; i < seqLen; i++) {
      embeddings.push(TP.getRow(weights.embedding, tokenIds[i]));
    }

    // Step 3: Positional encoding
    var positionalEncodings = computePositionalEncoding(seqLen, dModel);
    var positionedEmbeddings = TP.matAdd(embeddings, positionalEncodings);

    // Step 4 & 5: Multi-head attention
    var X = positionedEmbeddings;
    var qkvPerHead = [];
    var headOutputs = [];

    for (var h = 0; h < nHeads; h++) {
      var headResult = computeHeadAttention(X, weights.heads[h], seqLen, dK, causal);
      qkvPerHead.push(headResult);
      headOutputs.push(headResult.headOutput);
    }

    // Concatenate head outputs and multiply by W_O
    var concatenated = TP.concatCols(headOutputs);
    var multiHeadOutput = TP.matMul(concatenated, weights.W_O);

    // Step 6: Residual connection + LayerNorm (post-attention)
    var postResidual1 = TP.matAdd(X, multiHeadOutput);
    var postLayerNorm1 = TP.layerNormMatrix(postResidual1, weights.gamma1, weights.beta1);

    // Step 7: Feed-Forward Network
    // FFN(x) = ReLU(x * W1 + b1) * W2 + b2
    var ffnHidden = TP.relu(TP.addBias(TP.matMul(postLayerNorm1, weights.W1), weights.b1));
    var ffnOutput = TP.addBias(TP.matMul(ffnHidden, weights.W2), weights.b2);

    // Residual connection + LayerNorm (post-FFN)
    var postResidual2 = TP.matAdd(postLayerNorm1, ffnOutput);
    var postLayerNorm2 = TP.layerNormMatrix(postResidual2, weights.gamma2, weights.beta2);

    // Step 8: Output projection
    // For DNA (decoder): use last token's representation
    // For protein (encoder): use masked position's representation
    var outputIdx;
    if (mode === 'dna') {
      outputIdx = seqLen - 1; // last position predicts next token
    } else {
      outputIdx = Math.floor(seqLen / 2); // masked position
    }

    var outputRep = TP.getRow(postLayerNorm2, outputIdx);

    // logits = h * W_vocab + b_vocab
    // outputRep is 1D (length 8), W_vocab is 8 x vocabSize
    var vocabSize = weights.W_vocab[0].length;
    var logits = [];
    for (var v = 0; v < vocabSize; v++) {
      var logit = weights.b_vocab[v];
      for (var d = 0; d < dModel; d++) {
        logit += outputRep[d] * weights.W_vocab[d][v];
      }
      logits.push(logit);
    }

    // Softmax → probabilities
    var probabilities = TP.softmax(logits);

    // Argmax → predicted token
    var predictedIdx = TP.argmax(probabilities);
    var vocabRev = mode === 'dna' ? TP.DNA_VOCAB_REV : TP.PROTEIN_VOCAB_REV;
    var predictedToken = vocabRev[predictedIdx];

    return {
      tokenIds: tokenIds,
      embeddings: embeddings,
      positionalEncodings: positionalEncodings,
      positionedEmbeddings: positionedEmbeddings,
      qkvPerHead: qkvPerHead,
      multiHeadOutput: multiHeadOutput,
      postResidual1: postResidual1,
      postLayerNorm1: postLayerNorm1,
      ffnHidden: ffnHidden,
      ffnOutput: ffnOutput,
      postResidual2: postResidual2,
      postLayerNorm2: postLayerNorm2,
      outputIdx: outputIdx,
      logits: logits,
      probabilities: probabilities,
      predictedToken: predictedToken
    };
  };

  // Also expose the PE computation for visualization in Step 3
  TP.computePositionalEncoding = computePositionalEncoding;

})();
