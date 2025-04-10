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
            1835, 852, 262, 4494, 1933, 1323, 3301, 3348, 1800, 2575, 1369, 4363, 1799, 3030, 3178,
            3500, 1691, 4626, 589, 848, 2710, 1069, 2127, 3824, 2520, 3304, 758, 4197, 3306, 2489,
            3248, 2975, 1924, 3835, 3138, 1376, 905, 4323, 1665, 3888, 394, 2324, 1752, 3475, 2613,
            3782, 805, 5091, 3294, 2989, 1110, 2255, 3230, 4298, 2628, 3715, 1963, 2363, 322, 2861,
            3106, 5041, 2766, 3626, 1612, 3614, 1020, 1256, 2058, 3689, 2962, 3328, 1623, 5438, 42,
            3031, 2904, 2578, 1867, 3380, 1314, 2607, 1522, 5372, 185, 3568, 2305, 3470, 1631,
            1253, 251, 2210, 743, 2367, 2414, 3151, 2101, 1695, 314, 3150, 1911, 3446, 2518, 4622,
            2760, 5564, 394, 2605, 667, 4233, 2072, 967, 1087, 2017, 1321, 4448, 691, 5145, 568,
            2847, 195, 2282, 747, 4824, 1109, 624, 1343, 5250, 1131, 3543, 1691, 2270, 1268, 1767,
            676, 1678, 2526, 3797, 2042, 3418, 2380, 4049, 1794, 5700, 1735, 4583, 2565, 4815,
            1880, 5966, 2366, 5894, 3109, 3078, 1284, 3658, 2181, 2205, 1004, 3218, 2759, 5332, 21,
            4822, 2616, 1549, 627, 4894, 2589, 4510, 3226, 2310, 2182, 5329, 1710, 2256, 1185,
            2427, 436, 2935, 2377, 3550, 1851, 1700, 1296, 3981, 1251, 2173, 1452, 3058, 145, 3106,
            2738, 3109, 126, 3304, 2747, 1372, 2561, 2135, 2629, 3274, 755, 4183, 2423, 3907, 463,
            3490, 1100, 2122, 925, 3819, 2629, 4102, 2390, 1750, 1507, 869, 1921, 5432, 621, 3013,
            1202, 1860, 176, 1571, 144, 1886, 569, 3097, 1087, 3013, 2932, 2770, 2857, 2182, 1587,
            4112, 1853, 3722, 2843, 5523, 3004, 4087, 95, 1260, 2812, 3768, 688, 2139, 785, 3719,
            1169, 2507, 1538, 4459,
        ];
        assert_eq!(expect_fg_hat, fg_hat);
    }
}
