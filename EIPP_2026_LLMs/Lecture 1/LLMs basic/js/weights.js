// weights.js — Rigged weight matrices for DNA and protein modes
// Depends on math-utils.js constants. Stored in TP.WEIGHTS keyed by mode.

(function() {
  'use strict';
  var TP = window.TP;

  // Helper: create a seeded pseudo-random number generator for reproducibility
  function seededRandom(seed) {
    var s = seed;
    return function() {
      s = (s * 1664525 + 1013904223) & 0xffffffff;
      return (s >>> 0) / 0xffffffff;
    };
  }

  // Helper: create a matrix filled with small random values
  function randomMatrix(rows, cols, rng, scale) {
    scale = scale || 0.1;
    var mat = [];
    for (var i = 0; i < rows; i++) {
      var row = [];
      for (var j = 0; j < cols; j++) {
        row.push((rng() - 0.5) * 2 * scale);
      }
      mat.push(row);
    }
    return mat;
  }

  // Helper: create a random bias vector
  function randomBias(n, rng, scale) {
    scale = scale || 0.05;
    var b = [];
    for (var i = 0; i < n; i++) {
      b.push((rng() - 0.5) * 2 * scale);
    }
    return b;
  }

  // =========================================================================
  // DNA WEIGHTS
  // =========================================================================
  function buildDNAWeights() {
    var rng = seededRandom(42);
    var D = TP.D_MODEL; // 8
    var DK = TP.D_K;    // 2
    var DFF = TP.D_FF;  // 32
    var VOCAB = 4;       // A, T, G, C

    // --- Embedding matrix (4 x 8) ---
    // A and T are near-opposite (complementary base pair)
    // G and C are near-opposite (complementary base pair)
    // A and G share some similarity (purines), T and C share some (pyrimidines)
    var embedding = [
      // A (purine):  strong positive in dims 0-3, mild in 4-7
      [ 0.8,  0.6,  0.4,  0.3,  0.2,  0.1, -0.1, -0.2],
      // T (pyrimidine, complement of A): near-opposite of A
      [-0.7, -0.5, -0.3, -0.4,  0.3, -0.1,  0.2,  0.1],
      // G (purine):  similar direction to A in some dims, different in others
      [ 0.5,  0.3, -0.6, -0.4, -0.3,  0.7,  0.4,  0.2],
      // C (pyrimidine, complement of G): near-opposite of G
      [-0.4, -0.2,  0.5,  0.5,  0.4, -0.6, -0.3, -0.1]
    ];

    // --- Attention weights for 4 heads (each W_Q, W_K, W_V is 8 x 2) ---
    var heads = [];

    // Head 0: Complementary base pairing
    // Q and K projections designed so that A's query matches T's key (and vice versa)
    // and G's query matches C's key
    heads.push({
      W_Q: [
        [ 0.6,  0.1],  // dim 0
        [-0.5,  0.1],  // dim 1
        [ 0.1,  0.6],  // dim 2
        [ 0.1, -0.5],  // dim 3
        [ 0.0,  0.0],  // dim 4
        [ 0.0,  0.0],  // dim 5
        [ 0.0,  0.0],  // dim 6
        [ 0.0,  0.0]   // dim 7
      ],
      W_K: [
        [-0.5,  0.1],  // dim 0: flipped sign from Q → A's Q matches T's K
        [ 0.6,  0.1],  // dim 1
        [ 0.1, -0.5],  // dim 2: flipped → G's Q matches C's K
        [ 0.1,  0.6],  // dim 3
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.0,  0.0]
      ],
      W_V: randomMatrix(D, DK, rng, 0.3)
    });

    // Head 1: Local/sequential context (±1 position)
    // This is achieved through positional encoding interaction, but we boost it
    // by making Q and K focus on dimensions 4-7 where PE varies most locally
    heads.push({
      W_Q: [
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.7,  0.0],  // dim 4: high-frequency PE component
        [ 0.0,  0.7],  // dim 5
        [ 0.4,  0.0],  // dim 6
        [ 0.0,  0.4]   // dim 7
      ],
      W_K: [
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.7,  0.0],
        [ 0.0,  0.7],
        [ 0.4,  0.0],
        [ 0.0,  0.4]
      ],
      W_V: randomMatrix(D, DK, rng, 0.3)
    });

    // Head 2: Purine-pyrimidine patterns
    // Q/K designed so purines (A,G) attend to each other, pyrimidines (T,C) attend to each other
    heads.push({
      W_Q: [
        [ 0.5,  0.2],  // dims 0-1 carry purine signal
        [ 0.3,  0.4],
        [-0.3,  0.2],  // dims 2-3 carry pyrimidine signal (opposite sign)
        [-0.2,  0.3],
        [ 0.1,  0.0],
        [ 0.2,  0.1],
        [ 0.0,  0.1],
        [ 0.1,  0.0]
      ],
      W_K: [
        [ 0.5,  0.2],
        [ 0.3,  0.4],
        [-0.3,  0.2],
        [-0.2,  0.3],
        [ 0.1,  0.0],
        [ 0.2,  0.1],
        [ 0.0,  0.1],
        [ 0.1,  0.0]
      ],
      W_V: randomMatrix(D, DK, rng, 0.3)
    });

    // Head 3: Diffuse/distributed attention
    heads.push({
      W_Q: randomMatrix(D, DK, rng, 0.2),
      W_K: randomMatrix(D, DK, rng, 0.2),
      W_V: randomMatrix(D, DK, rng, 0.3)
    });

    // --- W_O: Output projection (8 x 8) ---
    // Maps concatenated 4-head outputs (4*2=8) back to d_model=8
    // Use near-identity with some mixing
    var W_O = [
      [ 0.8,  0.1,  0.0,  0.0,  0.1,  0.0,  0.0,  0.0],
      [ 0.1,  0.8,  0.0,  0.0,  0.0,  0.1,  0.0,  0.0],
      [ 0.0,  0.0,  0.8,  0.1,  0.0,  0.0,  0.1,  0.0],
      [ 0.0,  0.0,  0.1,  0.8,  0.0,  0.0,  0.0,  0.1],
      [ 0.1,  0.0,  0.0,  0.0,  0.8,  0.1,  0.0,  0.0],
      [ 0.0,  0.1,  0.0,  0.0,  0.1,  0.8,  0.0,  0.0],
      [ 0.0,  0.0,  0.1,  0.0,  0.0,  0.0,  0.8,  0.1],
      [ 0.0,  0.0,  0.0,  0.1,  0.0,  0.0,  0.1,  0.8]
    ];

    // --- FFN weights ---
    // W1 (8 x 32), b1 (32), W2 (32 x 8), b2 (8)
    var W1 = randomMatrix(D, DFF, rng, 0.3);
    var b1 = randomBias(DFF, rng, 0.1);
    var W2 = randomMatrix(DFF, D, rng, 0.3);
    var b2 = randomBias(D, rng, 0.05);

    // --- Output projection (8 x 4) and bias (4) ---
    // Biased to predict complement of the input pattern:
    // If last token embedding points toward A → bias toward T, etc.
    // Dims 0-1 are strong for A, so negative weight on A, positive on T
    // Output weights calibrated against actual final representations
    // to predict complement of last base: A→T, T→A, G→C, C→G
    var W_vocab = [
      [-0.15,  0.14, -0.08,  0.09],  // dim 0
      [-0.11,  0.09, -0.02,  0.04],  // dim 1
      [-0.07,  0.06,  0.18, -0.18],  // dim 2
      [-0.05,  0.03,  0.12, -0.09],  // dim 3
      [ 0.10, -0.05,  0.09, -0.15],  // dim 4
      [ 0.01,  0.00, -0.16,  0.15],  // dim 5
      [ 0.11, -0.15, -0.13,  0.16],  // dim 6
      [ 0.17, -0.12, -0.01, -0.03]   // dim 7
    ];
    var b_vocab = [0.0, 0.0, 0.0, 0.0];

    // --- LayerNorm parameters (gamma=1, beta=0 effectively) ---
    var gamma = [1, 1, 1, 1, 1, 1, 1, 1];
    var beta  = [0, 0, 0, 0, 0, 0, 0, 0];

    return {
      embedding: embedding,
      heads: heads,
      W_O: W_O,
      W1: W1, b1: b1,
      W2: W2, b2: b2,
      W_vocab: W_vocab,
      b_vocab: b_vocab,
      gamma1: gamma.slice(), beta1: beta.slice(),
      gamma2: gamma.slice(), beta2: beta.slice()
    };
  }

  // =========================================================================
  // PROTEIN WEIGHTS
  // =========================================================================
  function buildProteinWeights() {
    var rng = seededRandom(137);
    var D = TP.D_MODEL; // 8
    var DK = TP.D_K;    // 2
    var DFF = TP.D_FF;  // 32
    var VOCAB = 20;

    // --- Embedding matrix (20 x 8) ---
    // Organized by physicochemical properties:
    // Hydrophobic: A(0), F(4), I(7), L(9), M(10), V(17), W(18) — cluster together
    // Charged: D(2), E(3), K(8), R(14) — cluster together
    // Small/polar: C(1), G(5), N(11), Q(13), S(15), T(16) — cluster together
    // Aromatic: F(4), H(6), W(18), Y(19) — partial overlap with hydrophobic
    // P(12) is special (proline, rigid)

    // Embedding design:
    // dims 0-1: hydrophobicity axis (positive = hydrophobic)
    // dims 2-3: charge axis (positive = positive charge, negative = negative charge)
    // dims 4-5: size axis (positive = large)
    // dims 6-7: polarity / special features

    var embedding = [
      // A(0): small, mildly hydrophobic
      [ 0.4,  0.3, -0.1,  0.0, -0.5, -0.4,  0.1,  0.0],
      // C(1): small, polar, special (disulfide)
      [ 0.1,  0.0, -0.1,  0.1, -0.4, -0.3,  0.5,  0.4],
      // D(2): negative charge
      [-0.5, -0.4, -0.7, -0.5,  0.0,  0.1,  0.2,  0.1],
      // E(3): negative charge, larger
      [-0.4, -0.3, -0.6, -0.5,  0.3,  0.2,  0.1,  0.0],
      // F(4): hydrophobic, aromatic, large
      [ 0.7,  0.6,  0.0,  0.1,  0.5,  0.4,  0.2,  0.1],
      // G(5): tiny, flexible
      [ 0.0,  0.1, -0.1,  0.0, -0.7, -0.6, -0.2, -0.1],
      // H(6): aromatic, can be charged
      [ 0.1,  0.0,  0.3,  0.2,  0.2,  0.3,  0.4,  0.3],
      // I(7): hydrophobic, branched
      [ 0.7,  0.5,  0.0,  0.0,  0.3,  0.2, -0.2, -0.1],
      // K(8): positive charge, long
      [-0.4, -0.3,  0.7,  0.6,  0.3,  0.2,  0.1,  0.0],
      // L(9): hydrophobic
      [ 0.7,  0.6,  0.0,  0.0,  0.2,  0.3, -0.1,  0.0],
      // M(10): hydrophobic, contains sulfur
      [ 0.6,  0.5,  0.0,  0.1,  0.3,  0.2,  0.3,  0.2],
      // N(11): polar, amide
      [-0.2, -0.1, -0.1,  0.1,  0.0,  0.1,  0.4,  0.3],
      // P(12): special (rigid, proline)
      [ 0.2,  0.1,  0.0,  0.0, -0.2, -0.1, -0.5, -0.6],
      // Q(13): polar, larger amide
      [-0.1,  0.0, -0.1,  0.2,  0.2,  0.1,  0.3,  0.4],
      // R(14): positive charge, very long
      [-0.3, -0.2,  0.6,  0.7,  0.4,  0.5,  0.2,  0.1],
      // S(15): small, polar
      [-0.1,  0.0, -0.1,  0.0, -0.4, -0.3,  0.3,  0.2],
      // T(16): small, polar, branched
      [ 0.0,  0.1, -0.1,  0.0, -0.3, -0.2,  0.3,  0.1],
      // V(17): hydrophobic, branched, small
      [ 0.6,  0.4,  0.0,  0.0,  0.0, -0.1, -0.2, -0.1],
      // W(18): hydrophobic, aromatic, largest
      [ 0.6,  0.5,  0.0,  0.1,  0.7,  0.6,  0.3,  0.2],
      // Y(19): aromatic, has hydroxyl
      [ 0.4,  0.3,  0.0,  0.1,  0.4,  0.5,  0.4,  0.3]
    ];

    // --- Attention heads ---
    var heads = [];

    // Head 0: Hydrophobic clustering
    // Q and K focus on dims 0-1 (hydrophobicity axis) so hydrophobic residues attend to each other
    heads.push({
      W_Q: [
        [ 0.7,  0.0],  // dim 0: hydrophobicity
        [ 0.0,  0.7],  // dim 1: hydrophobicity
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.1,  0.0],
        [ 0.0,  0.1],
        [ 0.0,  0.0],
        [ 0.0,  0.0]
      ],
      W_K: [
        [ 0.7,  0.0],
        [ 0.0,  0.7],
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.1,  0.0],
        [ 0.0,  0.1],
        [ 0.0,  0.0],
        [ 0.0,  0.0]
      ],
      W_V: randomMatrix(D, DK, rng, 0.3)
    });

    // Head 1: Charge interactions
    // Q and K focus on dims 2-3 (charge axis) so charged residues attend to each other
    heads.push({
      W_Q: [
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.7,  0.0],  // dim 2: charge
        [ 0.0,  0.7],  // dim 3: charge
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.1,  0.0],
        [ 0.0,  0.1]
      ],
      W_K: [
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.7,  0.0],
        [ 0.0,  0.7],
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.1,  0.0],
        [ 0.0,  0.1]
      ],
      W_V: randomMatrix(D, DK, rng, 0.3)
    });

    // Head 2: Local backbone context (±2 positions)
    // Focus on high-frequency positional encoding dimensions (4-7)
    heads.push({
      W_Q: [
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.6,  0.0],
        [ 0.0,  0.6],
        [ 0.4,  0.0],
        [ 0.0,  0.4]
      ],
      W_K: [
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.0,  0.0],
        [ 0.6,  0.0],
        [ 0.0,  0.6],
        [ 0.4,  0.0],
        [ 0.0,  0.4]
      ],
      W_V: randomMatrix(D, DK, rng, 0.3)
    });

    // Head 3: Diffuse/distributed
    heads.push({
      W_Q: randomMatrix(D, DK, rng, 0.2),
      W_K: randomMatrix(D, DK, rng, 0.2),
      W_V: randomMatrix(D, DK, rng, 0.3)
    });

    // --- W_O: Output projection (8 x 8) ---
    var W_O = [
      [ 0.8,  0.1,  0.0,  0.0,  0.1,  0.0,  0.0,  0.0],
      [ 0.1,  0.8,  0.0,  0.0,  0.0,  0.1,  0.0,  0.0],
      [ 0.0,  0.0,  0.8,  0.1,  0.0,  0.0,  0.1,  0.0],
      [ 0.0,  0.0,  0.1,  0.8,  0.0,  0.0,  0.0,  0.1],
      [ 0.1,  0.0,  0.0,  0.0,  0.8,  0.1,  0.0,  0.0],
      [ 0.0,  0.1,  0.0,  0.0,  0.1,  0.8,  0.0,  0.0],
      [ 0.0,  0.0,  0.1,  0.0,  0.0,  0.0,  0.8,  0.1],
      [ 0.0,  0.0,  0.0,  0.1,  0.0,  0.0,  0.1,  0.8]
    ];

    // --- FFN weights ---
    var W1 = randomMatrix(D, DFF, rng, 0.3);
    var b1 = randomBias(DFF, rng, 0.1);
    var W2 = randomMatrix(DFF, D, rng, 0.3);
    var b2 = randomBias(D, rng, 0.05);

    // --- Output projection (8 x 20) ---
    // Biased toward predicting residues from the same physicochemical class
    // Dims 0-1 = hydrophobicity → boost hydrophobic residues
    // Dims 2-3 = charge → boost charged residues
    // Dims 4-5 = size
    // Dims 6-7 = polarity

    var W_vocab = [];
    // For each of the 8 embedding dimensions, create a row of 20 weights
    // that maps the representation back to vocab space
    // We construct this to roughly invert the embedding logic
    for (var d = 0; d < D; d++) {
      var row = [];
      for (var v = 0; v < VOCAB; v++) {
        // Use the embedding values as a guide: if embedding[v][d] is positive,
        // a positive activation in dim d should boost that vocab item
        row.push(embedding[v][d] * 0.5 + (rng() - 0.5) * 0.1);
      }
      W_vocab.push(row);
    }

    var b_vocab = randomBias(VOCAB, rng, 0.05);

    // --- LayerNorm parameters ---
    var gamma = [1, 1, 1, 1, 1, 1, 1, 1];
    var beta  = [0, 0, 0, 0, 0, 0, 0, 0];

    return {
      embedding: embedding,
      heads: heads,
      W_O: W_O,
      W1: W1, b1: b1,
      W2: W2, b2: b2,
      W_vocab: W_vocab,
      b_vocab: b_vocab,
      gamma1: gamma.slice(), beta1: beta.slice(),
      gamma2: gamma.slice(), beta2: beta.slice()
    };
  }

  // Build and store both weight sets
  TP.WEIGHTS = {
    dna: buildDNAWeights(),
    protein: buildProteinWeights()
  };

})();
