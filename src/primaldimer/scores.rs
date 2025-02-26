pub static NN_SCORES: [[[[Option<f64>; 4]; 4]; 4]; 4] = [
    [
        [
            [None, None, None, Some(0.69)],
            [None, None, None, Some(1.33)],
            [None, None, None, Some(0.74)],
            [Some(0.61), Some(0.88), Some(0.14), Some(-1.0)],
        ],
        [
            [None, None, Some(0.17), None],
            [None, None, Some(0.47), None],
            [None, None, Some(-0.52), None],
            [Some(0.77), Some(1.33), Some(-1.44), Some(0.64)],
        ],
        [
            [None, Some(0.43), None, None],
            [None, Some(0.79), None, None],
            [None, Some(0.11), None, None],
            [Some(0.02), Some(-1.28), Some(-0.13), Some(0.71)],
        ],
        [
            [Some(0.61), None, None, None],
            [Some(0.77), None, None, None],
            [Some(0.02), None, None, None],
            [Some(-0.88), Some(0.73), Some(0.07), Some(0.69)],
        ],
    ],
    [
        [
            [None, None, None, Some(0.92)],
            [None, None, None, Some(1.05)],
            [Some(0.43), Some(0.75), Some(0.03), Some(-1.45)],
            [None, None, None, Some(0.75)],
        ],
        [
            [None, None, Some(0.81), None],
            [None, None, Some(0.79), None],
            [Some(0.79), Some(0.7), Some(-1.84), Some(0.62)],
            [None, None, Some(0.98), None],
        ],
        [
            [None, Some(0.75), None, None],
            [None, Some(0.7), None, None],
            [Some(0.11), Some(-2.17), Some(-0.11), Some(-0.47)],
            [None, Some(0.4), None, None],
        ],
        [
            [Some(0.88), None, None, None],
            [Some(1.33), None, None, None],
            [Some(-1.28), Some(0.4), Some(-0.32), Some(-0.12)],
            [Some(0.73), None, None, None],
        ],
    ],
    [
        [
            [None, None, None, Some(0.42)],
            [Some(0.17), Some(0.81), Some(-0.25), Some(-1.3)],
            [None, None, None, Some(0.44)],
            [None, None, None, Some(0.34)],
        ],
        [
            [None, None, Some(-0.25), None],
            [Some(0.47), Some(0.79), Some(-2.24), Some(0.62)],
            [None, None, Some(-1.11), None],
            [None, None, Some(-0.59), None],
        ],
        [
            [None, Some(0.03), None, None],
            [Some(-0.52), Some(-1.84), Some(-1.11), Some(0.08)],
            [None, Some(-0.11), None, None],
            [None, Some(-0.32), None, None],
        ],
        [
            [Some(0.14), None, None, None],
            [Some(-1.44), Some(0.98), Some(-0.59), Some(0.45)],
            [Some(-0.13), None, None, None],
            [Some(0.07), None, None, None],
        ],
    ],
    [
        [
            [Some(0.69), Some(0.92), Some(0.42), Some(-0.58)],
            [None, None, None, Some(0.97)],
            [None, None, None, Some(0.43)],
            [None, None, None, Some(0.68)],
        ],
        [
            [Some(1.33), Some(1.05), Some(-1.3), Some(0.97)],
            [None, None, Some(0.62), None],
            [None, None, Some(0.08), None],
            [None, None, Some(0.45), None],
        ],
        [
            [Some(0.74), Some(-1.45), Some(0.44), Some(0.43)],
            [None, Some(0.62), None, None],
            [None, Some(-0.47), None, None],
            [None, Some(-0.12), None, None],
        ],
        [
            [Some(-1.0), Some(0.75), Some(0.34), Some(0.68)],
            [Some(0.64), None, None, None],
            [Some(0.71), None, None, None],
            [Some(0.69), None, None, None],
        ],
    ],
];

pub static SEQ1_OVERHANG_ARRAY: [[[Option<f64>; 4]; 4]; 4] = [
    [
        [None, None, None, None],
        [None, None, None, None],
        [None, None, None, None],
        [Some(-0.51), Some(-0.42), Some(-0.62), Some(-0.71)],
    ],
    [
        [None, None, None, None],
        [None, None, None, None],
        [Some(-0.96), Some(-0.52), Some(-0.72), Some(-0.58)],
        [None, None, None, None],
    ],
    [
        [None, None, None, None],
        [Some(-0.58), Some(-0.34), Some(-0.56), Some(-0.61)],
        [None, None, None, None],
        [None, None, None, None],
    ],
    [
        [Some(-0.51), Some(-0.02), Some(0.48), Some(-0.1)],
        [None, None, None, None],
        [None, None, None, None],
        [None, None, None, None],
    ],
];

pub static SEQ2_OVERHANG_ARRAY: [[[Option<f64>; 4]; 4]; 4] = [
    [
        [None, None, None, None],
        [None, None, None, None],
        [None, None, None, None],
        [Some(-0.48), Some(-0.19), Some(-0.5), Some(-0.29)],
    ],
    [
        [None, None, None, None],
        [None, None, None, None],
        [Some(-0.92), Some(-0.23), Some(-0.44), Some(-0.35)],
        [None, None, None, None],
    ],
    [
        [None, None, None, None],
        [Some(-0.82), Some(-0.31), Some(-0.01), Some(-0.52)],
        [None, None, None, None],
        [None, None, None, None],
    ],
    [
        [Some(-0.12), Some(0.28), Some(-0.01), Some(0.13)],
        [None, None, None, None],
        [None, None, None, None],
        [None, None, None, None],
    ],
];

// Given two encodes bases; ie MATCH_ARRAY[base][base]
// With return if these bases are a match (bool)
pub static MATCH_ARRAY: [[bool; 4]; 4] = [
    [false, false, false, true],
    [false, false, true, false],
    [false, true, false, false],
    [true, false, false, false],
];
