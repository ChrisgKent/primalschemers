pub fn nn_dg_scores(seq1: &[u8], seq2: &[u8]) -> Option<f64> {
    match (seq1, seq2) {
        (b"AA", b"AT") => Some(0.69),
        (b"AA", b"CT") => Some(1.33),
        (b"AA", b"GT") => Some(0.74),
        (b"AA", b"TA") => Some(0.61),
        (b"AA", b"TC") => Some(0.88),
        (b"AA", b"TG") => Some(0.14),
        (b"AA", b"TT") => Some(-1.0),
        (b"AC", b"AG") => Some(0.17),
        (b"AC", b"CG") => Some(0.47),
        (b"AC", b"GG") => Some(-0.52),
        (b"AC", b"TA") => Some(0.77),
        (b"AC", b"TC") => Some(1.33),
        (b"AC", b"TG") => Some(-1.44),
        (b"AC", b"TT") => Some(0.64),
        (b"AG", b"AC") => Some(0.43),
        (b"AG", b"CC") => Some(0.79),
        (b"AG", b"GC") => Some(0.11),
        (b"AG", b"TA") => Some(0.02),
        (b"AG", b"TC") => Some(-1.28),
        (b"AG", b"TG") => Some(-0.13),
        (b"AG", b"TT") => Some(0.71),
        (b"AT", b"AA") => Some(0.61),
        (b"AT", b"CA") => Some(0.77),
        (b"AT", b"GA") => Some(0.02),
        (b"AT", b"TA") => Some(-0.88),
        (b"AT", b"TC") => Some(0.73),
        (b"AT", b"TG") => Some(0.07),
        (b"AT", b"TT") => Some(0.69),
        (b"CA", b"AT") => Some(0.92),
        (b"CA", b"CT") => Some(1.05),
        (b"CA", b"GA") => Some(0.43),
        (b"CA", b"GC") => Some(0.75),
        (b"CA", b"GG") => Some(0.03),
        (b"CA", b"GT") => Some(-1.45),
        (b"CA", b"TT") => Some(0.75),
        (b"CC", b"AG") => Some(0.81),
        (b"CC", b"CG") => Some(0.79),
        (b"CC", b"GA") => Some(0.79),
        (b"CC", b"GC") => Some(0.7),
        (b"CC", b"GG") => Some(-1.84),
        (b"CC", b"GT") => Some(0.62),
        (b"CC", b"TG") => Some(0.98),
        (b"CG", b"AC") => Some(0.75),
        (b"CG", b"CC") => Some(0.7),
        (b"CG", b"GA") => Some(0.11),
        (b"CG", b"GC") => Some(-2.17),
        (b"CG", b"GG") => Some(-0.11),
        (b"CG", b"GT") => Some(-0.47),
        (b"CG", b"TC") => Some(0.4),
        (b"CT", b"AA") => Some(0.88),
        (b"CT", b"CA") => Some(1.33),
        (b"CT", b"GA") => Some(-1.28),
        (b"CT", b"GC") => Some(0.4),
        (b"CT", b"GG") => Some(-0.32),
        (b"CT", b"GT") => Some(-0.12),
        (b"CT", b"TA") => Some(0.73),
        (b"GA", b"AT") => Some(0.42),
        (b"GA", b"CA") => Some(0.17),
        (b"GA", b"CC") => Some(0.81),
        (b"GA", b"CG") => Some(-0.25),
        (b"GA", b"CT") => Some(-1.3),
        (b"GA", b"GT") => Some(0.44),
        (b"GA", b"TT") => Some(0.34),
        (b"GC", b"AG") => Some(-0.25),
        (b"GC", b"CA") => Some(0.47),
        (b"GC", b"CC") => Some(0.79),
        (b"GC", b"CG") => Some(-2.24),
        (b"GC", b"CT") => Some(0.62),
        (b"GC", b"GG") => Some(-1.11),
        (b"GC", b"TG") => Some(-0.59),
        (b"GG", b"AC") => Some(0.03),
        (b"GG", b"CA") => Some(-0.52),
        (b"GG", b"CC") => Some(-1.84),
        (b"GG", b"CG") => Some(-1.11),
        (b"GG", b"CT") => Some(0.08),
        (b"GG", b"GC") => Some(-0.11),
        (b"GG", b"TC") => Some(-0.32),
        (b"GT", b"AA") => Some(0.14),
        (b"GT", b"CA") => Some(-1.44),
        (b"GT", b"CC") => Some(0.98),
        (b"GT", b"CG") => Some(-0.59),
        (b"GT", b"CT") => Some(0.45),
        (b"GT", b"GA") => Some(-0.13),
        (b"GT", b"TA") => Some(0.07),
        (b"TA", b"AA") => Some(0.69),
        (b"TA", b"AC") => Some(0.92),
        (b"TA", b"AG") => Some(0.42),
        (b"TA", b"AT") => Some(-0.58),
        (b"TA", b"CT") => Some(0.97),
        (b"TA", b"GT") => Some(0.43),
        (b"TA", b"TT") => Some(0.68),
        (b"TC", b"AA") => Some(1.33),
        (b"TC", b"AC") => Some(1.05),
        (b"TC", b"AG") => Some(-1.3),
        (b"TC", b"AT") => Some(0.97),
        (b"TC", b"CG") => Some(0.62),
        (b"TC", b"GG") => Some(0.08),
        (b"TC", b"TG") => Some(0.45),
        (b"TG", b"AA") => Some(0.74),
        (b"TG", b"AC") => Some(-1.45),
        (b"TG", b"AG") => Some(0.44),
        (b"TG", b"AT") => Some(0.43),
        (b"TG", b"CC") => Some(0.62),
        (b"TG", b"GC") => Some(-0.47),
        (b"TG", b"TC") => Some(-0.12),
        (b"TT", b"AA") => Some(-1.0),
        (b"TT", b"AC") => Some(0.75),
        (b"TT", b"AG") => Some(0.34),
        (b"TT", b"AT") => Some(0.68),
        (b"TT", b"CA") => Some(0.64),
        (b"TT", b"GA") => Some(0.71),
        (b"TT", b"TA") => Some(0.69),
        _ => None,
    }
}

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

pub fn seq1_overhang_dg(seq1: &u8, seq2: &u8, seq1_overhang: &u8) -> Option<f64> {
    match (seq1, seq2, seq1_overhang) {
        (b'A', b'T', b'A') => Some(-0.51),
        (b'C', b'G', b'A') => Some(-0.96),
        (b'G', b'C', b'A') => Some(-0.58),
        (b'T', b'A', b'A') => Some(-0.51),
        (b'A', b'T', b'C') => Some(-0.42),
        (b'C', b'G', b'C') => Some(-0.52),
        (b'G', b'C', b'C') => Some(-0.34),
        (b'T', b'A', b'C') => Some(-0.02),
        (b'A', b'T', b'G') => Some(-0.62),
        (b'C', b'G', b'G') => Some(-0.72),
        (b'G', b'C', b'G') => Some(-0.56),
        (b'T', b'A', b'G') => Some(0.48),
        (b'A', b'T', b'T') => Some(-0.71),
        (b'C', b'G', b'T') => Some(-0.58),
        (b'G', b'C', b'T') => Some(-0.61),
        (b'T', b'A', b'T') => Some(-0.1),
        _ => None,
    }
}

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

pub fn seq2_overhang_dg(seq1: &u8, seq2: &u8, seq2_overhang: &u8) -> Option<f64> {
    match (seq1, seq2, seq2_overhang) {
        (b'A', b'T', b'A') => Some(-0.48),
        (b'C', b'G', b'A') => Some(-0.92),
        (b'G', b'C', b'A') => Some(-0.82),
        (b'T', b'A', b'A') => Some(-0.12),
        (b'A', b'T', b'C') => Some(-0.19),
        (b'C', b'G', b'C') => Some(-0.23),
        (b'G', b'C', b'C') => Some(-0.31),
        (b'T', b'A', b'C') => Some(0.28),
        (b'A', b'T', b'G') => Some(-0.50),
        (b'C', b'G', b'G') => Some(-0.44),
        (b'G', b'C', b'G') => Some(-0.01),
        (b'T', b'A', b'G') => Some(-0.01),
        (b'A', b'T', b'T') => Some(-0.29),
        (b'C', b'G', b'T') => Some(-0.35),
        (b'G', b'C', b'T') => Some(-0.52),
        (b'T', b'A', b'T') => Some(0.13),
        _ => None,
    }
}

// Given two encodes bases; ie MATCH_ARRAY[base][base]
// With return if these bases are a match (bool)
pub static MATCH_ARRAY: [[bool; 4]; 4] = [
    [false, false, false, true],
    [false, false, true, false],
    [false, true, false, false],
    [true, false, false, false],
];

pub fn match_array(base1: u8, base2: u8) -> bool {
    // Given two encodes bases; ie MATCH_ARRAY[base][base]
    // With return if these bases are a match (bool)
    match (base1, base2) {
        (b'A', b'T') => true,
        (b'T', b'A') => true,
        (b'C', b'G') => true,
        (b'G', b'C') => true,
        _ => false,
    }
}
