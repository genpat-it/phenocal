//! Reader for allelic profiles (`--profiles`): distances are computed here,
//! so a cgMLST/wgMLST allele-call table can be given directly instead of a
//! precomputed distance matrix.
//!
//! The table has one row per sample and one column per locus; the first
//! header cell is any label (`FILE` for chewBBACA `results_alleles.tsv`).
//! Calls are compared as in cgmlst-dists, chewBBACA and ReporTree: the
//! distance between two samples is the number of loci called in both whose
//! alleles differ, so missing calls never count as differences.
//!
//! A call is read as follows:
//! - a leading chewBBACA `INF-` is removed (`INF-42` is allele 42);
//! - an integer is an allele id, compared by value (`007` = `7`); `0` is missing;
//! - empty cells, `-`, `?`, `N`, `NA` and the chewBBACA classes (`LNF`,
//!   `PLOT3`, `PLOT5`, `LOTSC`, `NIPH`, `NIPHEM`, `ALM`, `ASM`, `PAMA`, `EXC`,
//!   `INV`) are missing;
//! - any other text (allele hashes, `*`-prefixed novel alleles) is an allele
//!   compared as text.

use crate::pairs::match_name;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

const MISSING: &[&str] = &[
    "", "-", "?", "N", "NA", "LNF", "PLOT3", "PLOT5", "LOTSC", "NIPH", "NIPHEM", "ALM", "ASM",
    "PAMA", "EXC", "INV",
];

/// Allele key of one call: `None` when missing.
fn allele_key(cell: &str) -> Option<AlleleKey<'_>> {
    let s = cell.trim();
    let s = s.strip_prefix("INF-").unwrap_or(s);
    if !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()) {
        let t = s.trim_start_matches('0');
        return if t.is_empty() {
            None
        } else {
            Some(AlleleKey::Int(t))
        };
    }
    if MISSING.iter().any(|m| m.eq_ignore_ascii_case(s)) {
        None
    } else {
        Some(AlleleKey::Text(s))
    }
}

/// Integers are keyed without leading zeros, so they compare by value at any
/// length; text is kept apart from integers.
#[derive(Hash, PartialEq, Eq)]
enum AlleleKey<'a> {
    Int(&'a str),
    Text(&'a str),
}

pub struct ProfileStats {
    pub n_samples: usize,
    pub n_loci: usize,
    pub n_matched: usize,
    pub n_matched_by_extension: usize,
    pub n_matched_by_prefix: usize,
    pub prefix_example: Option<(String, String)>,
    /// missing calls among the samples with a phenotype, out of n_matched * n_loci
    pub n_missing_calls: usize,
    pub n_pairs: usize,
}

/// Count of loci called in both profiles whose alleles differ (0 = missing).
fn hamming(a: &[u32], b: &[u32]) -> u32 {
    a.iter()
        .zip(b)
        .map(|(&x, &y)| ((x != y) & (x != 0) & (y != 0)) as u32)
        .sum()
}

/// Read an allelic-profile table, keep the samples present in `known` (the
/// phenotype table, sample -> cohort), and call `emit(sample_i, sample_j, d)`
/// for every pair of them from different cohorts, with `i` before `j` in file
/// order. Distances are computed on `threads` threads; the emitted order does
/// not depend on it.
pub fn read_profiles<F: FnMut(&str, &str, f64)>(
    path: &str,
    known: &HashMap<String, String>,
    threads: usize,
    mut emit: F,
) -> Result<ProfileStats, String> {
    let file = File::open(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let mut lines = BufReader::new(file).lines();
    let header = match lines.next() {
        Some(l) => l.map_err(|e| format!("cannot read {path}: {e}"))?,
        None => return Err(format!("{path}: empty file")),
    };
    let n_loci = header.trim_end_matches(['\r', '\n']).split('\t').count() - 1;
    if n_loci == 0 {
        return Err(format!(
            "{path}: header has no locus columns; expected a sample column followed \
             by one column per locus"
        ));
    }

    // allele text -> per-locus integer id (1, 2, ...; 0 is missing)
    let mut ids: Vec<HashMap<String, u32>> = vec![HashMap::new(); n_loci];
    let mut names: Vec<String> = Vec::new();
    let mut cohorts: Vec<&str> = Vec::new();
    let mut calls: Vec<u32> = Vec::new();
    let mut n_samples = 0usize;
    let (mut by_ext, mut by_prefix, mut n_missing) = (0usize, 0usize, 0usize);
    let mut prefix_example = None;
    let mut seen: HashMap<String, String> = HashMap::new();
    for line in lines {
        let line = line.map_err(|e| format!("cannot read {path}: {e}"))?;
        let line = line.trim_end_matches(['\r', '\n']);
        if line.trim().is_empty() {
            continue;
        }
        n_samples += 1;
        let mut cells = line.split('\t');
        let name = cells.next().unwrap_or("").trim();
        let id = match match_name(name, known).map_err(|e| format!("{path}: {e}"))? {
            Some((id, rule)) => {
                if rule == 1 {
                    by_ext += 1;
                } else if rule == 2 {
                    by_prefix += 1;
                    if prefix_example.is_none() {
                        prefix_example = Some((name.to_string(), id.clone()));
                    }
                }
                id
            }
            None => continue,
        };
        if let Some(prev) = seen.insert(id.clone(), name.to_string()) {
            return Err(format!(
                "{path}: rows '{prev}' and '{name}' both match sample '{id}'"
            ));
        }
        let start = calls.len();
        for (l, cell) in cells.enumerate() {
            if l == n_loci {
                return Err(format!(
                    "{path}: row '{name}' has more calls than the {n_loci} loci in the header"
                ));
            }
            let v = match allele_key(cell) {
                None => {
                    n_missing += 1;
                    0
                }
                Some(key) => {
                    let text = match key {
                        AlleleKey::Int(t) => t.to_string(),
                        AlleleKey::Text(t) => format!("\u{0}{t}"),
                    };
                    let next = ids[l].len() as u32 + 1;
                    *ids[l].entry(text).or_insert(next)
                }
            };
            calls.push(v);
        }
        if calls.len() - start != n_loci {
            return Err(format!(
                "{path}: row '{name}' has {} calls, expected {n_loci}",
                calls.len() - start
            ));
        }
        cohorts.push(known[&id].as_str());
        names.push(id);
    }
    let n = names.len();
    if n == 0 {
        let ex_p: Vec<&str> = known.keys().take(3).map(|s| s.as_str()).collect();
        return Err(format!(
            "{path}: no profile matches a phenotype sample (phenotypes: {ex_p:?}). \
             Profile names must equal the sample ids, possibly followed by '_...' or \
             a FASTA extension."
        ));
    }

    // cross-cohort distances, rows dealt round-robin to the threads so the
    // work is balanced over the triangle
    let threads = threads.clamp(1, n);
    let row = |i: usize| &calls[i * n_loci..(i + 1) * n_loci];
    let per_thread: Vec<Vec<(u32, u32, u32)>> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..threads)
            .map(|t| {
                let (row, cohorts) = (&row, &cohorts);
                s.spawn(move || {
                    let mut out = Vec::new();
                    for i in (t..n).step_by(threads) {
                        let a = row(i);
                        for j in i + 1..n {
                            if cohorts[i] != cohorts[j] {
                                out.push((i as u32, j as u32, hamming(a, row(j))));
                            }
                        }
                    }
                    out
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut by_row: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n];
    for (i, j, d) in per_thread.into_iter().flatten() {
        by_row[i as usize].push((j, d));
    }
    let mut n_pairs = 0usize;
    for (i, js) in by_row.iter().enumerate() {
        for &(j, d) in js {
            emit(&names[i], &names[j as usize], d as f64);
            n_pairs += 1;
        }
    }
    Ok(ProfileStats {
        n_samples,
        n_loci,
        n_matched: n,
        n_matched_by_extension: by_ext,
        n_matched_by_prefix: by_prefix,
        prefix_example,
        n_missing_calls: n_missing,
        n_pairs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(name: &str, text: &str) -> String {
        let p = std::env::temp_dir().join(format!("phenocal_prof_{}_{name}", std::process::id()));
        std::fs::write(&p, text).unwrap();
        p.to_string_lossy().into_owned()
    }

    fn known(pairs: &[(&str, &str)]) -> HashMap<String, String> {
        pairs
            .iter()
            .map(|(s, c)| (s.to_string(), c.to_string()))
            .collect()
    }

    fn collect(
        path: &str,
        k: &HashMap<String, String>,
        threads: usize,
    ) -> Result<Vec<(String, String, f64)>, String> {
        let mut v = Vec::new();
        read_profiles(path, k, threads, |a, b, d| {
            v.push((a.to_string(), b.to_string(), d))
        })?;
        Ok(v)
    }

    #[test]
    fn reads_calls_like_cgmlst_dists() {
        assert!(allele_key("0").is_none());
        assert!(allele_key("000").is_none());
        assert!(allele_key("LNF").is_none());
        assert!(allele_key("PLOT3").is_none());
        assert!(allele_key("-").is_none());
        assert!(allele_key("").is_none());
        assert!(allele_key("INF-42") == Some(AlleleKey::Int("42")));
        assert!(allele_key("007") == Some(AlleleKey::Int("7")));
        assert!(allele_key("*486") == Some(AlleleKey::Text("*486")));
        assert_eq!(hamming(&[0, 2, 3], &[1, 2, 4]), 1);
        assert_eq!(hamming(&[0, 0], &[1, 2]), 0);
    }

    #[test]
    fn counts_differences_on_loci_called_in_both() {
        let p = write(
            "basic.tsv",
            "FILE\tl1\tl2\tl3\tl4\n\
             a.fasta\t1\tINF-5\tLNF\t*9\n\
             b_shovill.fasta\t1\t5\t2\t*9\n\
             c\t2\t6\t3\t9\n\
             x\t1\t1\t1\t1\n",
        );
        let k = known(&[("a", "L1"), ("b", "L2"), ("c", "L2")]);
        let got = collect(&p, &k, 1).unwrap();
        // a-b: l1 same, l2 INF-5 = 5, l3 missing, l4 same -> 0
        // a-c: l1, l2 differ, l3 missing, l4 '*9' != '9' -> 3
        // b-c: same cohort -> not emitted; x has no phenotype
        assert_eq!(
            got,
            vec![("a".into(), "b".into(), 0.0), ("a".into(), "c".into(), 3.0)]
        );
        for t in 2..5 {
            assert_eq!(collect(&p, &k, t).unwrap(), got);
        }
    }

    #[test]
    fn rejects_ragged_rows_and_duplicate_samples() {
        let k = known(&[("a", "L1"), ("b", "L2")]);
        let p = write("ragged.tsv", "FILE\tl1\tl2\na\t1\t2\nb\t1\n");
        assert!(collect(&p, &k, 1).unwrap_err().contains("expected 2"));
        let p = write("dup.tsv", "FILE\tl1\na\t1\na.fasta\t2\nb\t1\n");
        assert!(collect(&p, &k, 1).unwrap_err().contains("both match"));
        let p = write("none.tsv", "FILE\tl1\nz\t1\n");
        assert!(collect(&p, &k, 1)
            .unwrap_err()
            .contains("no profile matches"));
    }
}
