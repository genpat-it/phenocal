//! Tiny dense linear algebra: solve a symmetric positive-(semi)definite system
//! by Gaussian elimination with partial pivoting. Used for the K x K weighted
//! least-squares normal equations (K = number of free cohorts), so it is small.

/// Solve `A x = b` for `x`, where `A` is `n x n` row-major. Returns `None` if
/// the system is singular (e.g. a disconnected calibration graph).
pub fn solve(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut m = a.to_vec();
    let mut y = b.to_vec();
    for col in 0..n {
        // partial pivot
        let mut piv = col;
        let mut best = m[col * n + col].abs();
        for r in (col + 1)..n {
            let v = m[r * n + col].abs();
            if v > best {
                best = v;
                piv = r;
            }
        }
        if best < 1e-12 {
            return None; // singular
        }
        if piv != col {
            for c in 0..n {
                m.swap(col * n + c, piv * n + c);
            }
            y.swap(col, piv);
        }
        let d = m[col * n + col];
        for r in (col + 1)..n {
            let f = m[r * n + col] / d;
            if f == 0.0 {
                continue;
            }
            for c in col..n {
                m[r * n + c] -= f * m[col * n + c];
            }
            y[r] -= f * y[col];
        }
    }
    // back-substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for c in (i + 1)..n {
            s -= m[i * n + c] * x[c];
        }
        x[i] = s / m[i * n + i];
    }
    Some(x)
}

/// Median of a slice (clones+sorts; inputs here are small).
pub fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 {
        s[n / 2]
    } else {
        0.5 * (s[n / 2 - 1] + s[n / 2])
    }
}

/// Quantile (linear interpolation) of a slice, q in [0,1].
pub fn quantile(sorted: &[f64], q: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    if n == 1 {
        return sorted[0];
    }
    let pos = q * (n as f64 - 1.0);
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    let frac = pos - lo as f64;
    sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solve_2x2_known() {
        // [[2,1],[1,3]] x = [3,5] -> x = [0.8, 1.4]
        let x = solve(&[2.0, 1.0, 1.0, 3.0], &[3.0, 5.0], 2).unwrap();
        assert!((x[0] - 0.8).abs() < 1e-9);
        assert!((x[1] - 1.4).abs() < 1e-9);
    }

    #[test]
    fn solve_singular_is_none() {
        // rank-1 matrix -> singular
        assert!(solve(&[1.0, 2.0, 2.0, 4.0], &[1.0, 2.0], 2).is_none());
    }

    #[test]
    fn median_odd_and_even() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[1.0, 2.0, 3.0, 4.0]), 2.5);
    }

    #[test]
    fn quantile_endpoints_and_interpolation() {
        let s = [0.0, 1.0, 2.0, 3.0, 4.0];
        assert_eq!(quantile(&s, 0.0), 0.0);
        assert_eq!(quantile(&s, 1.0), 4.0);
        assert_eq!(quantile(&s, 0.5), 2.0);
    }
}
