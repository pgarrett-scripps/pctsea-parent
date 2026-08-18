#[cfg(test)]
pub(crate) fn pearson(left: &[f64], right: &[f64]) -> Option<f64> {
    if left.len() != right.len() || left.len() < 2 {
        return None;
    }
    let n = left.len() as f64;
    let mean_left = left.iter().sum::<f64>() / n;
    let mean_right = right.iter().sum::<f64>() / n;
    let mut covariance = 0.0;
    let mut variance_left = 0.0;
    let mut variance_right = 0.0;
    for (&a, &b) in left.iter().zip(right) {
        let da = a - mean_left;
        let db = b - mean_right;
        covariance += da * db;
        variance_left += da * da;
        variance_right += db * db;
    }
    let denominator = (variance_left * variance_right).sqrt();
    (denominator > 0.0).then_some(covariance / denominator)
}

#[cfg(test)]
pub(crate) fn cosine(left: &[f64], right: &[f64]) -> Option<f64> {
    if left.len() != right.len() || left.is_empty() {
        return None;
    }
    let dot = left.iter().zip(right).map(|(a, b)| a * b).sum::<f64>();
    let norm_left = left.iter().map(|x| x * x).sum::<f64>().sqrt();
    let norm_right = right.iter().map(|x| x * x).sum::<f64>().sqrt();
    let denominator = norm_left * norm_right;
    (denominator > 0.0).then_some(dot / denominator)
}

#[cfg(test)]
pub(crate) fn dot(left: &[f64], right: &[f64]) -> Option<f64> {
    if left.len() != right.len() || left.is_empty() {
        return None;
    }
    Some(left.iter().zip(right).map(|(a, b)| a * b).sum())
}

pub(crate) fn hypergeometric_upper_tail(
    population: usize,
    successes: usize,
    draws: usize,
    observed: usize,
) -> f64 {
    if population == 0 || successes > population || draws > population {
        return f64::NAN;
    }
    let maximum = successes.min(draws);
    let minimum = draws.saturating_sub(population - successes);
    if observed <= minimum {
        return 1.0;
    }
    if observed > maximum {
        return 0.0;
    }
    let logs = log_factorials(population);
    let denominator = log_choose(&logs, population, draws);
    let terms: Vec<f64> = (observed..=maximum)
        .map(|k| {
            log_choose(&logs, successes, k) + log_choose(&logs, population - successes, draws - k)
                - denominator
        })
        .collect();
    let peak = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    (peak.exp() * terms.iter().map(|term| (term - peak).exp()).sum::<f64>()).min(1.0)
}

fn log_factorials(maximum: usize) -> Vec<f64> {
    let mut values = vec![0.0; maximum + 1];
    for value in 2..=maximum {
        values[value] = values[value - 1] + (value as f64).ln();
    }
    values
}

fn log_choose(logs: &[f64], n: usize, k: usize) -> f64 {
    if k > n {
        f64::NEG_INFINITY
    } else {
        logs[n] - logs[k] - logs[n - k]
    }
}

pub(crate) fn benjamini_hochberg(values: &[f64]) -> Vec<f64> {
    let mut indexed: Vec<_> = values
        .iter()
        .copied()
        .enumerate()
        .filter(|(_, value)| value.is_finite())
        .collect();
    indexed.sort_by(|left, right| left.1.total_cmp(&right.1));
    let count = indexed.len();
    let mut result = vec![f64::NAN; values.len()];
    let mut previous = 1.0_f64;
    for (rank, &(original, value)) in indexed.iter().enumerate().rev() {
        let adjusted = (value * count as f64 / (rank + 1) as f64).min(previous);
        previous = adjusted;
        result[original] = adjusted.min(1.0);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn correlations_are_sane() {
        assert!((pearson(&[1.0, 2.0, 3.0], &[2.0, 4.0, 6.0]).unwrap() - 1.0).abs() < 1e-12);
        assert!((cosine(&[1.0, 0.0], &[1.0, 0.0]).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn hypergeometric_known_case() {
        let p = hypergeometric_upper_tail(20, 7, 5, 3);
        assert!((p - 0.206_785_345_7).abs() < 1e-9, "{p}");
    }

    #[test]
    fn bh_preserves_order_and_monotonicity() {
        let adjusted = benjamini_hochberg(&[0.01, 0.04, 0.03, f64::NAN]);
        assert_eq!(&adjusted[..3], &[0.03, 0.04, 0.04]);
        assert!(adjusted[3].is_nan());
    }
}
