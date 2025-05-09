// These can be manifest constants, because there is a table of zetas which is computed from them, which is also constant.
pub const Q: i32 = 8380417;
pub const ZETA: i32 = 1753;
pub const ZETAS: [i32; 256] = [
    0, 4808194, 3765607, 3761513, 5178923, 5496691, 5234739, 5178987, 7778734, 3542485, 2682288,
    2129892, 3764867, 7375178, 557458, 7159240, 5010068, 4317364, 2663378, 6705802, 4855975,
    7946292, 676590, 7044481, 5152541, 1714295, 2453983, 1460718, 7737789, 4795319, 2815639,
    2283733, 3602218, 3182878, 2740543, 4793971, 5269599, 2101410, 3704823, 1159875, 394148,
    928749, 1095468, 4874037, 2071829, 4361428, 3241972, 2156050, 3415069, 1759347, 7562881,
    4805951, 3756790, 6444618, 6663429, 4430364, 5483103, 3192354, 556856, 3870317, 2917338,
    1853806, 3345963, 1858416, 3073009, 1277625, 5744944, 3852015, 4183372, 5157610, 5258977,
    8106357, 2508980, 2028118, 1937570, 4564692, 2811291, 5396636, 7270901, 4158088, 1528066,
    482649, 1148858, 5418153, 7814814, 169688, 2462444, 5046034, 4213992, 4892034, 1987814,
    5183169, 1736313, 235407, 5130263, 3258457, 5801164, 1787943, 5989328, 6125690, 3482206,
    4197502, 7080401, 6018354, 7062739, 2461387, 3035980, 621164, 3901472, 7153756, 2925816,
    3374250, 1356448, 5604662, 2683270, 5601629, 4912752, 2312838, 7727142, 7921254, 348812,
    8052569, 1011223, 6026202, 4561790, 6458164, 6143691, 1744507, 1753, 6444997, 5720892, 6924527,
    2660408, 6600190, 8321269, 2772600, 1182243, 87208, 636927, 4415111, 4423672, 6084020, 5095502,
    4663471, 8352605, 822541, 1009365, 5926272, 6400920, 1596822, 4423473, 4620952, 6695264,
    4969849, 2678278, 4611469, 4829411, 635956, 8129971, 5925040, 4234153, 6607829, 2192938,
    6653329, 2387513, 4768667, 8111961, 5199961, 3747250, 2296099, 1239911, 4541938, 3195676,
    2642980, 1254190, 8368000, 2998219, 141835, 8291116, 2513018, 7025525, 613238, 7070156,
    6161950, 7921677, 6458423, 4040196, 4908348, 2039144, 6500539, 7561656, 6201452, 6757063,
    2105286, 6006015, 6346610, 586241, 7200804, 527981, 5637006, 6903432, 1994046, 2491325,
    6987258, 507927, 7192532, 7655613, 6545891, 5346675, 8041997, 2647994, 3009748, 5767564,
    4148469, 749577, 4357667, 3980599, 2569011, 6764887, 1723229, 1665318, 2028038, 1163598,
    5011144, 3994671, 8368538, 7009900, 3020393, 3363542, 214880, 545376, 7609976, 3105558,
    7277073, 508145, 7826699, 860144, 3430436, 140244, 6866265, 6195333, 3123762, 2358373, 6187330,
    5365997, 6663603, 2926054, 7987710, 8077412, 3531229, 4405932, 4606686, 1900052, 7598542,
    1054478, 7648983,
];

// posimod
pub fn mod_q_pos(a: i32) -> i32 {
    let mut v = a;
    while v < 0 {
        v += Q;
    }
    v % Q
}

// mod+/-
pub fn mod_pm(a: i32, q: i32) -> i32 {
    let r = a % q;
    if r > q >> 1 {
        r - q
    } else if r < (-(q >> 1)) {
        r + q
    } else {
        r
    }
}

pub fn mod_vec(v: &[i32; 256], q: i32) -> [i32; 256] {
    let mut rv = [0; 256];
    for idx in 0..256 {
        rv[idx] = mod_pm(v[idx], q);
    }
    rv
}

pub fn power2_round(r: i32, d: i32) -> (i32, i32) {
    let r_plus = r % Q;
    let r_0 = mod_pm(r_plus, 1 << d);
    (((r_plus as i32) - r_0) >> d, r_0)
}

pub fn decompose(r: i32, gamma_2: i32) -> (i32, i32) {
    let r_plus = mod_q_pos(r);
    let mut r_0 = mod_pm(r_plus, gamma_2 << 1);
    let r_1 = if r_plus - r_0 == (Q - 1) {
        r_0 -= 1;
        0
    } else {
        (r_plus - r_0) / (gamma_2 << 1)
    };
    (r_1, r_0)
}

pub fn high_bits(r: i32, gamma_2: i32) -> i32 {
    decompose(r, gamma_2).0
}

pub fn high_bits_vec(r: &[i32; 256], gamma_2: i32) -> [i32; 256] {
    let mut result: [i32; 256] = [0; 256];
    for (i, v) in r.iter().enumerate() {
        result[i] = high_bits(*v, gamma_2)
    }
    result
}

pub fn low_bits(r: i32, gamma_2: i32) -> i32 {
    decompose(r, gamma_2).1
}

pub fn low_bits_vec(r: &[i32; 256], gamma_2: i32) -> [i32; 256] {
    let mut result: [i32; 256] = [0; 256];
    for (i, v) in r.iter().enumerate() {
        result[i] = low_bits(*v, gamma_2)
    }
    result
}

pub fn make_hint(z: i32, r: i32, gamma_2: i32) -> bool {
    let r_1 = high_bits(r, gamma_2);
    let v_1 = high_bits(r + z, gamma_2);
    r_1 != v_1
}

pub fn use_hint(h: bool, r: i32, gamma_2: i32) -> i32 {
    let m = (Q - 1) / (gamma_2 << 1);
    let (r_1, r_0) = decompose(r, gamma_2);
    if h && r_0 > 0 {
        (r_1 + 1) % (m as i32)
    } else if h && r_0 <= 0 {
        (r_1 - 1) % (m as i32)
    } else {
        r_1
    }
}

/// Multiply a and b and reduce the result mod Q. Could be more efficient.
/// @todo make it more efficient :-)
pub fn red_mul(a: i32, b: i32) -> i32 {
    (((a as i64) * (b as i64)) % (Q as i64)) as i32
}

pub fn poly_power2_round(a: &[i32; 256], d: i32) -> ([i32; 256], [i32; 256]) {
    let mut r_1: [i32; 256] = [0; 256];
    let mut r_2: [i32; 256] = [0; 256];
    a.iter().enumerate().for_each(|(idx, val)| {
        (r_1[idx], r_2[idx]) = power2_round(*val, d);
    });
    (r_1, r_2)
}

pub fn poly_add(a: &[i32; 256], b: &[i32; 256]) -> [i32; 256] {
    let mut rv = [0; 256];
    for i in 0..256 {
        rv[i] = (a[i] + b[i]) % Q;
    }
    rv
}

pub fn poly_sub(a: &[i32; 256], b: &[i32; 256]) -> [i32; 256] {
    let mut rv = [0; 256];
    for i in 0..256 {
        rv[i] = (a[i] - b[i]) % Q; // mod_pm(a[i] - b[i], Q);
    }
    rv
}

pub fn poly_mulntt(a: &[i32; 256], b: &[i32; 256]) -> [i32; 256] {
    let mut rv = [0; 256];
    for i in 0..256 {
        rv[i] = red_mul(a[i], b[i]);
    }
    rv
}

pub fn poly_count_ones(a: &[i32; 256]) -> usize {
    a.iter().filter(|x| **x == 1).count()
}

pub fn poly_make_hint(a: &[i32; 256], b: &[i32; 256], gamma_2: i32) -> [i32; 256] {
    let mut rv = [0; 256];
    for i in 0..256 {
        rv[i] = match make_hint(a[i], b[i], gamma_2) {
            true => 1,
            false => 0,
        };
    }
    rv
}

pub fn poly_pm(a: &[i32; 256]) -> [i32; 256] {
    let mut rv = [0; 256];
    for i in 0..256 {
        rv[i] = mod_pm(a[i], Q);
    }
    rv
}

pub fn minus(a: &[i32; 256]) -> [i32; 256] {
    let mut rv = [0; 256];
    for i in 0..256 {
        rv[i] = -a[i];
    }
    rv
}

#[cfg(test)]
mod tests {
    //use super::*;
    // @todo there are no vectors for the above, so let's leave them for now and we'll test them as part of the overall test
    // vectors and maybe put some regression tests here when that happens.
}
