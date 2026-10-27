//! Readers for the genomic-distance input (`--pairs`).
//!
//! Two layouts are accepted:
//! - **long**: one row per pair, header columns `sample_i  sample_j  distance`;
//! - **matrix**: a square distance matrix, as written by chewBBACA
//!   (`distance_matrix.tsv`, upper-triangular, or `distance_matrix_symmetric.tsv`),
//!   cgmlst-dists, cgDist or any tool with the same TSV layout. The first header cell is a label (`FILE`
//!   for chewBBACA); the other header cells, and the first cell of each row,
//!   are the sample names, rows in the same order as columns. cgMLST and
//!   wgMLST matrices have the same layout.
//!
//! The matrix may be symmetric, or triangular with the other half empty
//! (chewBBACA) or filled with zeros (cgmlst-dists `--matrix-format
//! lower-tri`/`upper-tri`); each unordered pair is read once.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub const LONG_I: &[&str] = &["sample_i", "acc_i", "i"];
pub const LONG_J: &[&str] = &["sample_j", "acc_j", "j"];
pub const LONG_D: &[&str] = &["distance", "dist", "d"];

#[derive(Debug, PartialEq)]
pub enum Layout {
    Long,
    Matrix,
}

/// Long layout when the header names the three long-format columns,
/// matrix layout otherwise.
pub fn detect(header: &str) -> Layout {
    let cols: Vec<String> = header
        .split('\t')
        .map(|h| h.trim().to_lowercase())
        .collect();
    let has = |aliases: &[&str]| cols.iter().any(|c| aliases.contains(&c.as_str()));
    if has(LONG_I) && has(LONG_J) && has(LONG_D) {
        Layout::Long
    } else {
        Layout::Matrix
    }
}

/// Remove a FASTA file extension: chewBBACA names genomes after their files.
pub fn strip_fasta_ext(name: &str) -> &str {
    let base = name.strip_suffix(".gz").unwrap_or(name);
    for ext in [".fasta", ".fna", ".fa", ".fas", ".fsa", ".ffn", ".fsta"] {
        if let Some(s) = base.strip_suffix(ext) {
            return s;
        }
    }
    name
}

/// Match a matrix sample name to a phenotype sample id: the name itself, the
/// name without its FASTA extension (`ERR1.fasta` -> `ERR1`), or the unique
/// phenotype id followed in the name by `_` or `.` (`ERR1_shovill.fasta` ->
/// `ERR1`, the assembler suffix chewBBACA keeps from the file name). Returns
/// the id and the rule used (0 exact, 1 extension, 2 prefix); an ambiguous
/// prefix is an error.
pub fn match_name<T>(
    name: &str,
    known: &HashMap<String, T>,
) -> Result<Option<(String, u8)>, String> {
    if known.contains_key(name) {
        return Ok(Some((name.to_string(), 0)));
    }
    let base = strip_fasta_ext(name);
    if base != name && known.contains_key(base) {
        return Ok(Some((base.to_string(), 1)));
    }
    let hits: Vec<&str> = base
        .char_indices()
        .filter(|&(i, c)| i > 0 && (c == '_' || c == '.'))
        .map(|(i, _)| &base[..i])
        .filter(|p| known.contains_key(*p))
        .collect();
    match hits.len() {
        0 => Ok(None),
        1 => Ok(Some((hits[0].to_string(), 2))),
        _ => Err(format!(
            "matrix sample '{name}' matches several phenotype samples {hits:?}; \
             rename it to the exact sample id"
        )),
    }
}

/// Which half of the matrix holds the distances.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Triangle {
    /// both halves, mirroring each other
    Symmetric,
    /// upper half; the lower one is empty (chewBBACA) or all zeros (cgmlst-dists)
    Upper,
    /// lower half; the upper one is empty or all zeros
    Lower,
}

pub struct MatrixStats {
    pub n_genomes: usize,
    pub n_matched: usize,
    pub n_matched_by_extension: usize,
    pub n_matched_by_prefix: usize,
    pub prefix_example: Option<(String, String)>,
    pub triangle: Triangle,
    pub n_missing_cells: usize,
}

fn is_missing(cell: &str) -> bool {
    cell.is_empty()
        || matches!(
            cell.to_lowercase().as_str(),
            "na" | "nan" | "null" | "none" | "-" | "."
        )
}

fn parse_cell(cell: &str, row: &str, col: &str) -> Result<Option<f64>, String> {
    let cell = cell.trim();
    if is_missing(cell) {
        return Ok(None);
    }
    cell.parse::<f64>()
        .map(Some)
        .map_err(|_| format!("non-numeric distance '{cell}' between '{row}' and '{col}'"))
}

fn read_header(path: &str) -> Result<Vec<String>, String> {
    let file = File::open(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let mut header = String::new();
    BufReader::new(file)
        .read_line(&mut header)
        .map_err(|e| format!("cannot read {path}: {e}"))?;
    let names: Vec<String> = header
        .trim_end_matches(['\r', '\n'])
        .split('\t')
        .skip(1)
        .map(|s| s.trim().to_string())
        .collect();
    if names.len() < 2 {
        return Err(format!(
            "{path}: header has {} sample column(s); expected either the long layout \
             (sample_i, sample_j, distance) or a square distance matrix",
            names.len()
        ));
    }
    Ok(names)
}

/// Call `f(row_index, cells)` for every data row, after checking that the
/// matrix is square and lists the samples in the same order in rows and columns.
fn for_each_row<F: FnMut(usize, &[&str]) -> Result<(), String>>(
    path: &str,
    names: &[String],
    mut f: F,
) -> Result<(), String> {
    let n = names.len();
    let file = File::open(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let mut r = 0usize;
    for line in BufReader::new(file).lines().skip(1) {
        let line = line.map_err(|e| format!("cannot read {path}: {e}"))?;
        let line = line.trim_end_matches(['\r', '\n']);
        if line.trim().is_empty() {
            continue;
        }
        if r == n {
            return Err(format!(
                "{path}: more rows than the {n} samples in the header"
            ));
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let row_name = cells[0].trim();
        if row_name != names[r] {
            return Err(format!(
                "{path}: row {} is '{row_name}' but column {} is '{}'; rows and columns \
                 must list the samples in the same order",
                r + 1,
                r + 1,
                names[r]
            ));
        }
        if cells.len() != n + 1 {
            return Err(format!(
                "{path}: row '{row_name}' has {} distance cells, expected {n}",
                cells.len() - 1
            ));
        }
        f(r, &cells[1..])?;
        r += 1;
    }
    if r != n {
        return Err(format!(
            "{path}: {r} rows for {n} samples in the header; the matrix must be square"
        ));
    }
    Ok(())
}

/// Read a square distance matrix and call `emit(sample_i, sample_j, d)` once
/// per unordered pair of samples present in `known` (the phenotype table).
/// Matrix names are matched to `known` with [`match_name`].
///
/// Two passes over the file, O(n) memory besides what `emit` keeps. The first
/// finds which half holds the distances: a symmetric matrix (each sample's
/// distances sum the same in both halves), or one half with the other empty
/// or all zeros, the two triangular layouts of chewBBACA and cgmlst-dists.
/// Any other matrix is rejected, so a zero-filled half is never read as
/// distances.
pub fn read_matrix<T, F: FnMut(&str, &str, f64)>(
    path: &str,
    known: &HashMap<String, T>,
    mut emit: F,
) -> Result<MatrixStats, String> {
    let names = read_header(path)?;
    let n = names.len();

    // matrix name -> phenotype sample id
    let mut resolved: Vec<Option<String>> = Vec::with_capacity(n);
    let mut by_ext = 0usize;
    let mut by_prefix = 0usize;
    let mut prefix_example = None;
    let mut seen: HashMap<String, usize> = HashMap::new();
    for (k, name) in names.iter().enumerate() {
        let id = match match_name(name, known).map_err(|e| format!("{path}: {e}"))? {
            Some((id, rule)) => {
                if rule == 1 {
                    by_ext += 1;
                } else if rule == 2 {
                    by_prefix += 1;
                    if prefix_example.is_none() {
                        prefix_example = Some((name.clone(), id.clone()));
                    }
                }
                Some(id)
            }
            None => None,
        };
        if let Some(id) = &id {
            if let Some(prev) = seen.insert(id.clone(), k) {
                return Err(format!(
                    "{path}: matrix columns '{}' and '{name}' both match sample '{id}'",
                    names[prev]
                ));
            }
        }
        resolved.push(id);
    }
    let n_matched = resolved.iter().filter(|r| r.is_some()).count();
    if n_matched == 0 {
        let ex_m: Vec<&str> = names.iter().take(3).map(|s| s.as_str()).collect();
        let ex_p: Vec<&str> = known.keys().take(3).map(|s| s.as_str()).collect();
        return Err(format!(
            "{path}: no matrix sample matches a phenotype sample (matrix: {ex_m:?}; \
             phenotypes: {ex_p:?}). Matrix names must equal the sample ids, possibly \
             followed by '_...' or a FASTA extension."
        ));
    }

    // pass 1: per-sample sums of each half, and whether a half has any value
    // or any non-zero value
    let mut sum_up = vec![0.0f64; n];
    let mut sum_lo = vec![0.0f64; n];
    let (mut up_any, mut up_nonzero, mut lo_any, mut lo_nonzero) = (false, false, false, false);
    for_each_row(path, &names, |r, cells| {
        for (c, cell) in cells.iter().enumerate() {
            if c == r {
                continue;
            }
            if let Some(d) = parse_cell(cell, &names[r], &names[c])? {
                let (sum, any, nonzero) = if c > r {
                    (&mut sum_up, &mut up_any, &mut up_nonzero)
                } else {
                    (&mut sum_lo, &mut lo_any, &mut lo_nonzero)
                };
                sum[r] += d;
                sum[c] += d;
                *any = true;
                *nonzero |= d != 0.0;
            }
        }
        Ok(())
    })?;
    let mirrored = (0..n).all(|k| {
        (sum_up[k] - sum_lo[k]).abs() <= 1e-9 * sum_up[k].abs().max(sum_lo[k].abs()).max(1.0)
    });
    let triangle = if !lo_any || (!lo_nonzero && up_nonzero) {
        Triangle::Upper
    } else if !up_any || (!up_nonzero && lo_nonzero) {
        Triangle::Lower
    } else if mirrored {
        Triangle::Symmetric
    } else {
        return Err(format!(
            "{path}: the two halves of the matrix differ, so it is neither symmetric \
             nor triangular"
        ));
    };

    // pass 2: emit each pair once, from the half that holds the distances
    let mut n_missing = 0usize;
    for_each_row(path, &names, |r, cells| {
        let range = match triangle {
            Triangle::Upper | Triangle::Symmetric => (r + 1)..n,
            Triangle::Lower => 0..r,
        };
        for c in range {
            match parse_cell(cells[c], &names[r], &names[c])? {
                Some(d) => {
                    let (i, j) = if c > r { (r, c) } else { (c, r) };
                    if let (Some(a), Some(b)) = (&resolved[i], &resolved[j]) {
                        emit(a, b, d);
                    }
                }
                None => n_missing += 1,
            }
        }
        Ok(())
    })?;
    Ok(MatrixStats {
        n_genomes: n,
        n_matched,
        n_matched_by_extension: by_ext,
        n_matched_by_prefix: by_prefix,
        prefix_example,
        triangle,
        n_missing_cells: n_missing,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(name: &str, text: &str) -> String {
        let p = std::env::temp_dir().join(format!("phenocal_pairs_{name}.tsv"));
        std::fs::write(&p, text).unwrap();
        p.to_str().unwrap().to_string()
    }

    fn known(ids: &[&str]) -> HashMap<String, ()> {
        ids.iter().map(|s| (s.to_string(), ())).collect()
    }

    fn collect(path: &str, k: &HashMap<String, ()>) -> Result<Vec<(String, String, f64)>, String> {
        let mut v = Vec::new();
        read_matrix(path, k, |a, b, d| v.push((a.into(), b.into(), d)))?;
        Ok(v)
    }

    #[test]
    fn detects_layout() {
        assert_eq!(detect("sample_i\tsample_j\tdistance"), Layout::Long);
        assert_eq!(detect("FILE\ta.fasta\tb.fasta"), Layout::Matrix);
        assert_eq!(detect("cgmlst-dists\ta\tb"), Layout::Matrix);
    }

    #[test]
    fn strips_fasta_extensions_only() {
        assert_eq!(strip_fasta_ext("ERR1_shovill.fasta"), "ERR1_shovill");
        assert_eq!(strip_fasta_ext("s1.fna.gz"), "s1");
        assert_eq!(strip_fasta_ext("s1.tsv.gz"), "s1.tsv.gz");
        assert_eq!(strip_fasta_ext("s1"), "s1");
    }

    // chewBBACA distance_matrix.tsv (upper-triangular) and
    // distance_matrix_symmetric.tsv for the same profiles: same pairs.
    #[test]
    fn triangular_and_symmetric_give_the_same_pairs() {
        let tri = write(
            "tri",
            "FILE\ta.fasta\tb.fasta\tc.fasta\n\
             a.fasta\t0\t3\t7\n\
             b.fasta\t\t0\t5\n\
             c.fasta\t\t\t0\n",
        );
        let sym = write(
            "sym",
            "FILE\ta.fasta\tb.fasta\tc.fasta\n\
             a.fasta\t0\t3\t7\n\
             b.fasta\t3\t0\t5\n\
             c.fasta\t7\t5\t0\n",
        );
        let k = known(&["a", "b", "c"]);
        let expect = vec![
            ("a".to_string(), "b".to_string(), 3.0),
            ("a".to_string(), "c".to_string(), 7.0),
            ("b".to_string(), "c".to_string(), 5.0),
        ];
        assert_eq!(collect(&tri, &k).unwrap(), expect);
        assert_eq!(collect(&sym, &k).unwrap(), expect);
    }

    #[test]
    fn matches_chewbbaca_file_names() {
        let k = known(&["ERR1", "ERR10", "s1_A"]);
        let m = |n: &str| match_name(n, &k).unwrap();
        assert_eq!(m("ERR1"), Some(("ERR1".into(), 0)));
        assert_eq!(m("ERR10.fasta"), Some(("ERR10".into(), 1)));
        assert_eq!(m("ERR10_shovill.fasta"), Some(("ERR10".into(), 2)));
        assert_eq!(m("ERR1_shovill"), Some(("ERR1".into(), 2)));
        assert_eq!(m("s1_A_spades.fna"), Some(("s1_A".into(), 2)));
        assert_eq!(m("ERR100_shovill.fasta"), None);
        let amb = known(&["s1", "s1_A"]);
        assert!(match_name("s1_A_x.fasta", &amb).is_err());
    }

    // cgmlst-dists --matrix-format lower-tri / upper-tri fill the other half with 0
    #[test]
    fn zero_filled_triangles_are_not_read_as_distances() {
        let k = known(&["a", "b", "c"]);
        let expect = vec![
            ("a".to_string(), "b".to_string(), 3.0),
            ("a".to_string(), "c".to_string(), 7.0),
            ("b".to_string(), "c".to_string(), 5.0),
        ];
        let lower = write(
            "lower0",
            "cgmlst-dists\ta\tb\tc\na\t0\t0\t0\nb\t3\t0\t0\nc\t7\t5\t0\n",
        );
        let upper = write(
            "upper0",
            "cgmlst-dists\ta\tb\tc\na\t0\t3\t7\nb\t0\t0\t5\nc\t0\t0\t0\n",
        );
        let mut got = collect(&lower, &k).unwrap();
        got.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert_eq!(got, expect);
        assert_eq!(collect(&upper, &k).unwrap(), expect);
    }

    #[test]
    fn skips_samples_without_phenotype() {
        let p = write(
            "partial",
            "FILE\ta\tb\tc\na\t0\t3\t7\nb\t\t0\t5\nc\t\t\t0\n",
        );
        let v = collect(&p, &known(&["a", "c"])).unwrap();
        assert_eq!(v, vec![("a".to_string(), "c".to_string(), 7.0)]);
    }

    #[test]
    fn rejects_asymmetric_and_misordered_matrices() {
        let asym = write(
            "asym",
            "FILE\ta\tb\tc\na\t0\t3\t7\nb\t4\t0\t5\nc\t7\t5\t0\n",
        );
        assert!(collect(&asym, &known(&["a", "b", "c"]))
            .unwrap_err()
            .contains("neither symmetric"));
        let order = write("order", "FILE\ta\tb\na\t0\t3\nc\t\t0\n");
        assert!(collect(&order, &known(&["a", "b"]))
            .unwrap_err()
            .contains("same order"));
        let none = write("none", "FILE\tx\ty\nx\t0\t1\ny\t\t0\n");
        assert!(collect(&none, &known(&["a", "b"]))
            .unwrap_err()
            .contains("no matrix sample matches"));
    }
}
