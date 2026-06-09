//! Minimal deterministic RNG (SplitMix64) + standard-normal sampler (Box-Muller).
//! Seeded for reproducible bootstrap draws without external crates.

pub struct Rng {
    state: u64,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng {
            state: seed.wrapping_add(0x9E3779B97F4A7C15),
        }
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }

    /// Uniform in (0, 1).
    #[inline]
    fn next_f64(&mut self) -> f64 {
        // 53-bit mantissa, shifted off zero
        let u = self.next_u64() >> 11;
        (u as f64 + 0.5) * (1.0 / 9007199254740992.0)
    }

    /// Uniform integer in [0, n) for bootstrap resampling.
    #[inline]
    pub fn next_index(&mut self, n: usize) -> usize {
        ((self.next_f64() * n as f64) as usize).min(n.saturating_sub(1))
    }

    /// One standard-normal draw via Box-Muller.
    #[inline]
    pub fn next_normal(&mut self) -> f64 {
        let u1 = self.next_f64();
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// One Student-t draw with `nu` degrees of freedom, rescaled to unit variance
    /// (so it is a drop-in heavier-tailed replacement for next_normal).
    /// t = Z / sqrt(W/nu), W ~ chi2(nu); Var(t)=nu/(nu-2) for nu>2, so multiply
    /// by sqrt((nu-2)/nu) to keep unit variance.
    #[inline]
    pub fn next_t(&mut self, nu: f64) -> f64 {
        let z = self.next_normal();
        let k = nu.round().max(1.0) as usize;
        let mut w = 0.0;
        for _ in 0..k {
            let g = self.next_normal();
            w += g * g;
        }
        let t = z / (w / nu).sqrt();
        if nu > 2.0 {
            t * ((nu - 2.0) / nu).sqrt()
        } else {
            t
        }
    }
}
