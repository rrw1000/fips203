use crate::basics;
use anyhow::Result;

fn base_case_multiply(a: (u32, u32), b: (u32, u32), gamma: u32) -> (u32, u32) {
    // this is all mod Q, so we need to reduce every two multiplies,
    let c0 = (((a.0 * b.0) % basics::Q) + ((a.1 * b.1 % basics::Q) * gamma)) % basics::Q;
    let c1 = (a.0 * b.1) % basics::Q + (a.1 * b.0) % basics::Q;
    (c0, c1)
}

pub fn multiply_ntts(f_bar: [u32; 256], g_bar: [u32; 256]) -> Result<[u32; 256]> {
    let mut h_bar: [u32; 256] = [0; 256];
    // We step through all the arrays, two entries at a time.
    // @todo could probably be cleverer about coding the degree-2 polys here rather than using the representation in the standard.
    for i in 0..128 {
        let zeta = basics::zeta_2(u32::try_from(i)?);
        // Probably as cheap to shift here as count by two and shift back, and clearer this way too.
        let i_2: usize = i << 1;
        (h_bar[i_2], h_bar[i_2 + 1]) = base_case_multiply(
            (f_bar[i_2], f_bar[i_2 + 1]),
            (g_bar[i_2], g_bar[i_2 + 1]),
            zeta,
        );
    }
    Ok(h_bar)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    #[test]
    fn test_base_case_multiply() {
        let (c0, c1) = base_case_multiply((325, 678), (914, 123), 293);
        // Calculated via the python REPL
        assert_eq!(c0, 351);
        assert_eq!(c1, 525);
    }

    #[test]
    fn test_multiply_ntts() {
        // This is "just" a regression test.
        let mut rng = rand::rngs::StdRng::seed_from_u64(925);
        let f_hat: [u32; 256] = (0..256)
            .map(|_| rng.random_range(0..basics::Q))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let g_hat: [u32; 256] = (0..256)
            .map(|_| rng.random_range(0..basics::Q))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let fg_hat = multiply_ntts(f_hat, g_hat).unwrap();
        let expect_fg_hat: [u32; 256] = [
            1835, 852, 1905, 4494, 3156, 1323, 1947, 3348, 1482, 2575, 3108, 4363, 1779, 3030,
            2121, 3500, 1907, 4626, 2687, 848, 2610, 1069, 2907, 3824, 2605, 3304, 419, 4197, 2513,
            2489, 3281, 2975, 2979, 3835, 2962, 1376, 1134, 4323, 1141, 3888, 1949, 2324, 362,
            3475, 428, 3782, 472, 5091, 2367, 2989, 218, 2255, 2549, 4298, 1826, 3715, 455, 2363,
            692, 2861, 1938, 5041, 3317, 3626, 3276, 3614, 3155, 1256, 1951, 3689, 3055, 3328,
            2886, 5438, 1928, 3031, 688, 2578, 1066, 3380, 3157, 2607, 2542, 5372, 459, 3568, 1228,
            3470, 910, 1253, 259, 2210, 1678, 2367, 2408, 3151, 2359, 1695, 1063, 3150, 2378, 3446,
            3294, 4622, 1824, 5564, 784, 2605, 3256, 4233, 1940, 967, 2021, 2017, 1339, 4448, 1758,
            5145, 684, 2847, 3185, 2282, 2421, 4824, 2278, 624, 1491, 5250, 2897, 3543, 1273, 2270,
            1974, 1767, 2876, 1678, 1140, 3797, 346, 3418, 502, 4049, 449, 5700, 1301, 4583, 3040,
            4815, 254, 5966, 2343, 5894, 3228, 3078, 2835, 3658, 266, 2205, 2779, 3218, 78, 5332,
            2614, 4822, 548, 1549, 2387, 4894, 2742, 4510, 1221, 2310, 796, 5329, 563, 2256, 1312,
            2427, 2005, 2935, 1111, 3550, 1359, 1700, 2197, 3981, 2754, 2173, 1179, 3058, 2946,
            3106, 1592, 3109, 1513, 3304, 1250, 1372, 2000, 2135, 1155, 3274, 2933, 4183, 1137,
            3907, 1647, 3490, 2735, 2122, 2385, 3819, 2865, 4102, 279, 1750, 26, 869, 3166, 5432,
            448, 3013, 505, 1860, 1156, 1571, 2408, 1886, 2009, 3097, 906, 3013, 515, 2770, 820,
            2182, 1091, 4112, 3096, 3722, 2385, 5523, 1761, 4087, 1189, 1260, 256, 3768, 2305,
            2139, 1811, 3719, 26, 2507, 3105, 4459,
        ];
        assert_eq!(expect_fg_hat, fg_hat);
    }
}
