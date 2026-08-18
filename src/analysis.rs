use crate::stats::{benjamini_hochberg, hypergeometric_upper_tail};
use crate::{AtlasData, Error, GeneQuery, LoadedQuery, Result};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::{Display, Formatter};
use std::str::FromStr;

/// Similarity used to rank single cells against the experimental query.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ScoringMethod {
    #[default]
    Pearson,
    Cosine,
    DotProduct,
}

impl Display for ScoringMethod {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Pearson => f.write_str("pearson"),
            Self::Cosine => f.write_str("cosine"),
            Self::DotProduct => f.write_str("dot-product"),
        }
    }
}

impl FromStr for ScoringMethod {
    type Err = Error;

    fn from_str(value: &str) -> Result<Self> {
        match value.trim().to_ascii_lowercase().as_str() {
            "pearson" | "correlation" => Ok(Self::Pearson),
            "cosine" | "normalized-dot-product" => Ok(Self::Cosine),
            "dot" | "dot-product" | "dot_product" => Ok(Self::DotProduct),
            _ => Err(Error::InvalidConfig(format!(
                "unknown scoring method {value:?}; use pearson, cosine, or dot-product"
            ))),
        }
    }
}

/// Controls scoring, filtering, and enrichment significance.
#[derive(Clone, Debug)]
pub struct AnalysisConfig {
    pub method: ScoringMethod,
    pub min_genes_per_cell: usize,
    /// Positive thresholds retain scores `>= threshold`; negative thresholds
    /// retain scores `<= threshold`. `None` retains every scorable cell.
    pub score_threshold: Option<f64>,
    pub permutations: usize,
    pub random_seed: u64,
    /// Empty means all datasets.
    pub datasets: Vec<String>,
    /// Match the original PCTSEA behavior by scoring only pairs where both
    /// experimental and single-cell values are positive.
    pub positive_pairs_only: bool,
}

impl Default for AnalysisConfig {
    fn default() -> Self {
        Self {
            method: ScoringMethod::Pearson,
            min_genes_per_cell: 3,
            score_threshold: Some(0.0),
            permutations: 1_000,
            random_seed: 0x05ee_d5ee_dd15_ca11,
            datasets: Vec::new(),
            positive_pairs_only: true,
        }
    }
}

/// A ranked single-cell score.
#[derive(Clone, Debug)]
pub struct CellScore {
    pub cell_index: usize,
    pub cell_name: String,
    pub cell_type: String,
    pub dataset: String,
    pub score: f64,
    pub genes_used: usize,
    pub passes_threshold: bool,
}

/// Enrichment statistics for one cell type.
#[derive(Clone, Debug)]
pub struct CellTypeResult {
    pub cell_type: String,
    pub cells_total: usize,
    pub cells_passing: usize,
    pub hypergeometric_p_value: f64,
    pub hypergeometric_fdr: f64,
    pub enrichment_score: f64,
    pub normalized_enrichment_score: f64,
    pub permutation_p_value: f64,
    pub permutation_fdr: f64,
    /// Genes most frequently present in passing cells of this type.
    pub top_genes: Vec<(String, usize)>,
}

/// Complete in-memory result of an analysis.
#[derive(Clone, Debug)]
pub struct AnalysisResult {
    pub matched_genes: Vec<String>,
    pub missing_genes: Vec<String>,
    pub cells_considered: usize,
    pub cells_scored: usize,
    pub cells_passing: usize,
    pub cell_scores: Vec<CellScore>,
    pub cell_types: Vec<CellTypeResult>,
}

/// Scores cells and tests cell-type enrichment against a file-backed atlas.
pub fn analyze<A: AtlasData + ?Sized>(
    atlas: &A,
    query: &GeneQuery,
    config: &AnalysisConfig,
) -> Result<AnalysisResult> {
    validate(config, query)?;
    let loaded = atlas.load_query(query)?;
    if config.min_genes_per_cell > loaded.gene_names.len() {
        return Err(Error::InvalidConfig(format!(
            "min_genes_per_cell ({}) exceeds the number of query genes found in the atlas ({})",
            config.min_genes_per_cell,
            loaded.gene_names.len()
        )));
    }
    let selected_datasets: HashSet<_> = config.datasets.iter().map(String::as_str).collect();
    let mut scores = score_cells(atlas, &loaded, config, &selected_datasets);
    let cells_considered = atlas
        .cells()
        .iter()
        .filter(|cell| {
            selected_datasets.is_empty() || selected_datasets.contains(cell.dataset.as_str())
        })
        .count();
    let cells_scored = scores.len();
    if cells_scored == 0 {
        return Err(Error::InvalidConfig(
            "no cells were scorable after dataset and minimum-gene filtering".into(),
        ));
    }
    scores.sort_by(|left, right| compare_scores(left.score, right.score, config.score_threshold));
    let passing_indices: Vec<_> = scores
        .iter()
        .enumerate()
        .filter_map(|(index, score)| score.passes_threshold.then_some(index))
        .collect();
    let cells_passing = passing_indices.len();

    let mut type_names: Vec<_> = scores
        .iter()
        .map(|score| score.cell_type.clone())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
    type_names.sort();
    let type_lookup: HashMap<_, _> = type_names
        .iter()
        .enumerate()
        .map(|(index, name)| (name.as_str(), index))
        .collect();
    let score_type_ids: Vec<_> = scores
        .iter()
        .map(|score| type_lookup[score.cell_type.as_str()])
        .collect();
    let ranked_types: Vec<_> = passing_indices
        .iter()
        .map(|&index| score_type_ids[index])
        .collect();
    let ranked_weights: Vec<_> = passing_indices
        .iter()
        .map(|&index| scores[index].score)
        .collect();

    let mut totals = vec![0_usize; type_names.len()];
    let mut passing = vec![0_usize; type_names.len()];
    for &type_id in &score_type_ids {
        totals[type_id] += 1;
    }
    for &type_id in &ranked_types {
        passing[type_id] += 1;
    }
    let real_enrichment = enrichment_scores(&ranked_types, &ranked_weights, type_names.len());
    let random_enrichment = permuted_enrichment(
        &ranked_types,
        &ranked_weights,
        type_names.len(),
        config.permutations,
        config.random_seed,
    );

    let hyper_p: Vec<_> = (0..type_names.len())
        .map(|type_id| {
            hypergeometric_upper_tail(
                cells_scored,
                totals[type_id],
                cells_passing,
                passing[type_id],
            )
        })
        .collect();
    let hyper_fdr = benjamini_hochberg(&hyper_p);
    let permutation_p: Vec<_> = real_enrichment
        .iter()
        .enumerate()
        .map(|(type_id, &real)| empirical_p_value(real, &random_enrichment[type_id]))
        .collect();
    let permutation_fdr = benjamini_hochberg(&permutation_p);
    let normalized: Vec<_> = real_enrichment
        .iter()
        .enumerate()
        .map(|(type_id, &real)| normalize_enrichment(real, &random_enrichment[type_id]))
        .collect();
    let top_genes = top_genes_by_type(
        &scores,
        &passing_indices,
        &score_type_ids,
        &loaded,
        type_names.len(),
    );

    let mut cell_types: Vec<_> = type_names
        .into_iter()
        .enumerate()
        .map(|(type_id, cell_type)| CellTypeResult {
            cell_type,
            cells_total: totals[type_id],
            cells_passing: passing[type_id],
            hypergeometric_p_value: hyper_p[type_id],
            hypergeometric_fdr: hyper_fdr[type_id],
            enrichment_score: real_enrichment[type_id],
            normalized_enrichment_score: normalized[type_id],
            permutation_p_value: permutation_p[type_id],
            permutation_fdr: permutation_fdr[type_id],
            top_genes: top_genes[type_id].clone(),
        })
        .collect();
    cell_types.sort_by(|left, right| {
        finite_first(left.permutation_fdr, right.permutation_fdr)
            .then_with(|| {
                right
                    .normalized_enrichment_score
                    .total_cmp(&left.normalized_enrichment_score)
            })
            .then_with(|| left.cell_type.cmp(&right.cell_type))
    });

    Ok(AnalysisResult {
        matched_genes: loaded.gene_names,
        missing_genes: loaded.missing_genes,
        cells_considered,
        cells_scored,
        cells_passing,
        cell_scores: scores,
        cell_types,
    })
}

fn validate(config: &AnalysisConfig, query: &GeneQuery) -> Result<()> {
    if config.min_genes_per_cell == 0 {
        return Err(Error::InvalidConfig(
            "min_genes_per_cell must be at least 1".into(),
        ));
    }
    if config.method == ScoringMethod::Pearson && config.min_genes_per_cell < 2 {
        return Err(Error::InvalidConfig(
            "Pearson scoring requires min_genes_per_cell of at least 2".into(),
        ));
    }
    if config.min_genes_per_cell > query.genes.len() {
        return Err(Error::InvalidConfig(format!(
            "min_genes_per_cell ({}) exceeds query size ({})",
            config.min_genes_per_cell,
            query.genes.len()
        )));
    }
    if config
        .score_threshold
        .is_some_and(|value| !value.is_finite())
    {
        return Err(Error::InvalidConfig(
            "score threshold must be finite".into(),
        ));
    }
    Ok(())
}

fn score_cells<A: AtlasData + ?Sized>(
    atlas: &A,
    loaded: &LoadedQuery,
    config: &AnalysisConfig,
    datasets: &HashSet<&str>,
) -> Vec<CellScore> {
    atlas
        .cells()
        .iter()
        .enumerate()
        .filter(|(_, cell)| datasets.is_empty() || datasets.contains(cell.dataset.as_str()))
        .filter_map(|(cell_index, cell)| {
            let (score, genes_used) = score_cell(
                &loaded.query_values,
                &loaded.cell_values[cell_index],
                config,
            )?;
            let passes_threshold = passes(score, config.score_threshold);
            Some(CellScore {
                cell_index,
                cell_name: cell.name.clone(),
                cell_type: cell.cell_type.clone(),
                dataset: cell.dataset.clone(),
                score,
                genes_used,
                passes_threshold,
            })
        })
        .collect()
}

fn score_cell(
    query_values: &[f64],
    expressions: &[f32],
    config: &AnalysisConfig,
) -> Option<(f64, usize)> {
    if query_values.len() != expressions.len() {
        return None;
    }
    let genes_used = included_pairs(query_values, expressions, config.positive_pairs_only).count();
    if genes_used < config.min_genes_per_cell {
        return None;
    }
    let score = match config.method {
        ScoringMethod::Pearson => {
            let count = genes_used as f64;
            let mean_query = included_pairs(query_values, expressions, config.positive_pairs_only)
                .map(|(query, _)| query)
                .sum::<f64>()
                / count;
            let mean_expression =
                included_pairs(query_values, expressions, config.positive_pairs_only)
                    .map(|(_, expression)| expression)
                    .sum::<f64>()
                    / count;
            let (covariance, query_variance, expression_variance) =
                included_pairs(query_values, expressions, config.positive_pairs_only).fold(
                    (0.0, 0.0, 0.0),
                    |(covariance, query_variance, expression_variance), (query, expression)| {
                        let query_delta = query - mean_query;
                        let expression_delta = expression - mean_expression;
                        (
                            covariance + query_delta * expression_delta,
                            query_variance + query_delta * query_delta,
                            expression_variance + expression_delta * expression_delta,
                        )
                    },
                );
            let denominator = (query_variance * expression_variance).sqrt();
            (denominator > 0.0).then_some(covariance / denominator)?
        }
        ScoringMethod::Cosine => {
            let dot = included_pairs(query_values, expressions, config.positive_pairs_only)
                .map(|(query, expression)| query * expression)
                .sum::<f64>();
            let query_norm = included_pairs(query_values, expressions, config.positive_pairs_only)
                .map(|(query, _)| query * query)
                .sum::<f64>()
                .sqrt();
            let expression_norm =
                included_pairs(query_values, expressions, config.positive_pairs_only)
                    .map(|(_, expression)| expression * expression)
                    .sum::<f64>()
                    .sqrt();
            let denominator = query_norm * expression_norm;
            (denominator > 0.0).then_some(dot / denominator)?
        }
        ScoringMethod::DotProduct => {
            included_pairs(query_values, expressions, config.positive_pairs_only)
                .map(|(query, expression)| query * expression)
                .sum()
        }
    };
    Some((score, genes_used))
}

fn included_pairs<'a>(
    query_values: &'a [f64],
    expressions: &'a [f32],
    positive_pairs_only: bool,
) -> impl Iterator<Item = (f64, f64)> + 'a {
    query_values
        .iter()
        .copied()
        .zip(expressions.iter().map(|&value| f64::from(value)))
        .filter(move |&(query, expression)| {
            !positive_pairs_only || (query > 0.0 && expression > 0.0)
        })
}

fn passes(score: f64, threshold: Option<f64>) -> bool {
    match threshold {
        None => true,
        Some(value) if value >= 0.0 => score >= value,
        Some(value) => score <= value,
    }
}

fn compare_scores(left: f64, right: f64, threshold: Option<f64>) -> std::cmp::Ordering {
    if threshold.is_some_and(|value| value < 0.0) {
        left.total_cmp(&right)
    } else {
        right.total_cmp(&left)
    }
}

/// Conventional weighted GSEA running-sum statistic over ranked passing cells.
#[cfg(test)]
fn enrichment_score(types: &[usize], scores: &[f64], target: usize) -> f64 {
    let hits = types.iter().filter(|&&type_id| type_id == target).count();
    if hits == 0 || hits == types.len() {
        return f64::NAN;
    }
    let weight_sum: f64 = types
        .iter()
        .zip(scores)
        .filter(|(type_id, _)| **type_id == target)
        .map(|(_, score)| score.abs())
        .sum();
    let weighted = weight_sum > 0.0;
    let miss_step = 1.0 / (types.len() - hits) as f64;
    let mut running = 0.0;
    let mut supremum = 0.0_f64;
    for (&type_id, &score) in types.iter().zip(scores) {
        if type_id == target {
            running += if weighted {
                score.abs() / weight_sum
            } else {
                1.0 / hits as f64
            };
        } else {
            running -= miss_step;
        }
        if running.abs() > supremum.abs() {
            supremum = running;
        }
    }
    supremum
}

/// Calculates the same running-sum enrichment statistic for every type in one
/// pass. For one target type, the running sum increases only at its hits and
/// decreases monotonically between hits. Its largest absolute value therefore
/// occurs immediately before or after a hit, or at the end of the ranking.
fn enrichment_scores(types: &[usize], scores: &[f64], type_count: usize) -> Vec<f64> {
    if types.len() != scores.len() || types.iter().any(|&type_id| type_id >= type_count) {
        return vec![f64::NAN; type_count];
    }
    let ranking_length = types.len();
    let mut hits = vec![0_usize; type_count];
    let mut weight_sums = vec![0.0; type_count];
    for (&type_id, &score) in types.iter().zip(scores) {
        hits[type_id] += 1;
        weight_sums[type_id] += score.abs();
    }
    let mut hits_seen = vec![0_usize; type_count];
    let mut hit_weight_seen = vec![0.0; type_count];
    let mut miss_steps = vec![f64::NAN; type_count];
    let mut supremums = vec![0.0_f64; type_count];
    for type_id in 0..type_count {
        if hits[type_id] == 0 || hits[type_id] == ranking_length {
            supremums[type_id] = f64::NAN;
        } else {
            miss_steps[type_id] = 1.0 / (ranking_length - hits[type_id]) as f64;
        }
    }
    for (position, (&type_id, &score)) in types.iter().zip(scores).enumerate() {
        if !supremums[type_id].is_finite() {
            continue;
        }
        let misses_before = position - hits_seen[type_id];
        let before_hit = hit_weight_seen[type_id] - misses_before as f64 * miss_steps[type_id];
        update_supremum(&mut supremums[type_id], before_hit);
        let hit_step = if weight_sums[type_id] > 0.0 {
            score.abs() / weight_sums[type_id]
        } else {
            1.0 / hits[type_id] as f64
        };
        hit_weight_seen[type_id] += hit_step;
        hits_seen[type_id] += 1;
        let after_hit = hit_weight_seen[type_id] - misses_before as f64 * miss_steps[type_id];
        update_supremum(&mut supremums[type_id], after_hit);
    }
    for type_id in 0..type_count {
        if !supremums[type_id].is_finite() {
            continue;
        }
        let misses = ranking_length - hits[type_id];
        let final_value = hit_weight_seen[type_id] - misses as f64 * miss_steps[type_id];
        update_supremum(&mut supremums[type_id], final_value);
    }
    supremums
}

fn update_supremum(supremum: &mut f64, candidate: f64) {
    // Mathematically equal positive and negative extrema can differ by one
    // floating-point rounding step. Retain the first extreme in ranking order.
    let scale = candidate.abs().max(supremum.abs()).max(1.0);
    let tie_tolerance = 16.0 * f64::EPSILON * scale;
    if candidate.abs() - supremum.abs() > tie_tolerance {
        *supremum = candidate;
    }
}

fn permuted_enrichment(
    types: &[usize],
    scores: &[f64],
    type_count: usize,
    permutations: usize,
    seed: u64,
) -> Vec<Vec<f64>> {
    let mut result = vec![Vec::with_capacity(permutations); type_count];
    if permutations == 0 || types.is_empty() {
        return result;
    }
    let mut rng = Rng::new(seed);
    let mut shuffled = types.to_vec();
    for _ in 0..permutations {
        shuffled.copy_from_slice(types);
        rng.shuffle(&mut shuffled);
        let enrichment = enrichment_scores(&shuffled, scores, type_count);
        for (values, score) in result.iter_mut().zip(enrichment) {
            values.push(score);
        }
    }
    result
}

fn empirical_p_value(real: f64, random: &[f64]) -> f64 {
    if !real.is_finite() || random.is_empty() {
        return f64::NAN;
    }
    let comparable: Vec<_> = random
        .iter()
        .copied()
        .filter(|value| value.is_finite() && value.signum() == real.signum())
        .collect();
    if comparable.is_empty() {
        return f64::NAN;
    }
    let extreme = comparable
        .iter()
        .filter(|&&value| {
            if real >= 0.0 {
                value >= real
            } else {
                value <= real
            }
        })
        .count();
    (extreme + 1) as f64 / (comparable.len() + 1) as f64
}

fn normalize_enrichment(real: f64, random: &[f64]) -> f64 {
    if !real.is_finite() {
        return f64::NAN;
    }
    let same_sign: Vec<_> = random
        .iter()
        .copied()
        .filter(|value| value.is_finite() && value.signum() == real.signum())
        .collect();
    if same_sign.is_empty() {
        return f64::NAN;
    }
    let expected = same_sign.iter().map(|value| value.abs()).sum::<f64>() / same_sign.len() as f64;
    if expected == 0.0 {
        f64::NAN
    } else {
        real / expected
    }
}

fn top_genes_by_type(
    scores: &[CellScore],
    passing_indices: &[usize],
    score_type_ids: &[usize],
    loaded: &LoadedQuery,
    type_count: usize,
) -> Vec<Vec<(String, usize)>> {
    let mut counts = vec![BTreeMap::<String, usize>::new(); type_count];
    for &score_index in passing_indices {
        let score = &scores[score_index];
        let type_id = score_type_ids[score_index];
        for (gene_index, &value) in loaded.cell_values[score.cell_index].iter().enumerate() {
            if value > 0.0 && loaded.query_values[gene_index] > 0.0 {
                *counts[type_id]
                    .entry(loaded.gene_names[gene_index].clone())
                    .or_default() += 1;
            }
        }
    }
    counts
        .into_iter()
        .map(|counts| {
            let mut values: Vec<_> = counts.into_iter().collect();
            values.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
            values.truncate(20);
            values
        })
        .collect()
}

fn finite_first(left: f64, right: f64) -> std::cmp::Ordering {
    match (left.is_finite(), right.is_finite()) {
        (true, true) => left.total_cmp(&right),
        (true, false) => std::cmp::Ordering::Less,
        (false, true) => std::cmp::Ordering::Greater,
        (false, false) => std::cmp::Ordering::Equal,
    }
}

struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(if seed == 0 {
            0x9e37_79b9_7f4a_7c15
        } else {
            seed
        })
    }

    fn next_u64(&mut self) -> u64 {
        let mut value = self.0;
        value ^= value >> 12;
        value ^= value << 25;
        value ^= value >> 27;
        self.0 = value;
        value.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    fn shuffle<T>(&mut self, values: &mut [T]) {
        for index in (1..values.len()).rev() {
            let other = (self.next_u64() % (index as u64 + 1)) as usize;
            values.swap(index, other);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::{cosine, dot, pearson};

    #[test]
    fn enrichment_detects_top_heavy_type() {
        let types = [0, 0, 1, 1];
        let scores = [4.0, 3.0, 2.0, 1.0];
        assert!(enrichment_score(&types, &scores, 0) > 0.9);
        assert!(enrichment_score(&types, &scores, 1) < -0.9);
    }

    fn assert_enrichment_equivalent(actual: f64, expected: f64) {
        if expected.is_nan() {
            assert!(actual.is_nan());
        } else if (actual - expected).abs() >= 1e-12 {
            assert!(
                (actual.abs() - expected.abs()).abs() < 1e-12,
                "expected={expected} actual={actual}"
            );
        }
    }

    #[test]
    fn permutations_are_deterministic() {
        let types = [0, 0, 1, 1, 1];
        let scores = [5.0, 4.0, 3.0, 2.0, 1.0];
        assert_eq!(
            permuted_enrichment(&types, &scores, 2, 20, 42),
            permuted_enrichment(&types, &scores, 2, 20, 42)
        );
    }

    #[test]
    fn single_pass_enrichment_matches_per_type_calculation() {
        let types = [0, 1, 0, 2, 1, 0, 2, 2, 1];
        let scores = [9.0, 8.0, 0.0, -6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
        let combined = enrichment_scores(&types, &scores, 4);
        for (type_id, &actual) in combined.iter().enumerate() {
            let expected = enrichment_score(&types, &scores, type_id);
            assert_enrichment_equivalent(actual, expected);
        }
    }

    #[test]
    fn single_pass_enrichment_matches_randomized_rankings() {
        let mut rng = Rng::new(91);
        for ranking_length in 2..80 {
            for type_count in 1..10 {
                let types: Vec<_> = (0..ranking_length)
                    .map(|_| (rng.next_u64() % type_count as u64) as usize)
                    .collect();
                let scores: Vec<_> = (0..ranking_length)
                    .map(|_| (rng.next_u64() % 10_000) as f64 / 100.0 - 50.0)
                    .collect();
                let combined = enrichment_scores(&types, &scores, type_count);
                for (type_id, &actual) in combined.iter().enumerate() {
                    let expected = enrichment_score(&types, &scores, type_id);
                    assert_enrichment_equivalent(actual, expected);
                }
            }
        }
    }

    #[test]
    fn allocation_free_cell_scores_match_slice_statistics() {
        let query = [1.0, 2.0, 3.0, 4.0];
        let expressions = [2.0_f32, 0.0, 4.0, 8.0];
        let filtered_query = [1.0, 3.0, 4.0];
        let filtered_expressions = [2.0, 4.0, 8.0];
        for (method, expected) in [
            (
                ScoringMethod::Pearson,
                pearson(&filtered_query, &filtered_expressions).unwrap(),
            ),
            (
                ScoringMethod::Cosine,
                cosine(&filtered_query, &filtered_expressions).unwrap(),
            ),
            (
                ScoringMethod::DotProduct,
                dot(&filtered_query, &filtered_expressions).unwrap(),
            ),
        ] {
            let config = AnalysisConfig {
                method,
                min_genes_per_cell: 3,
                ..AnalysisConfig::default()
            };
            let (actual, genes_used) = score_cell(&query, &expressions, &config).unwrap();
            assert_eq!(actual, expected);
            assert_eq!(genes_used, 3);
        }
    }
}
