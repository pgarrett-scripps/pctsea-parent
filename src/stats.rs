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

/// Returns average ranks and the sum of `t^3 - t` over tied groups.
pub(crate) fn average_ranks(values: &[f64]) -> Option<(Vec<f64>, f64)> {
    if values.is_empty() || values.iter().any(|value| !value.is_finite()) {
        return None;
    }
    let mut order: Vec<_> = (0..values.len()).collect();
    order.sort_by(|&left, &right| values[left].total_cmp(&values[right]));
    let mut ranks = vec![0.0; values.len()];
    let mut tie_sum = 0.0;
    let mut start = 0;
    while start < order.len() {
        let mut end = start + 1;
        while end < order.len() && values[order[end]] == values[order[start]] {
            end += 1;
        }
        let average = (start + 1 + end) as f64 / 2.0;
        for &index in &order[start..end] {
            ranks[index] = average;
        }
        let tied = (end - start) as f64;
        tie_sum += tied * tied * tied - tied;
        start = end;
    }
    Some((ranks, tie_sum))
}

pub(crate) fn standard_normal_two_sided_p(z: f64) -> f64 {
    if !z.is_finite() {
        return f64::NAN;
    }
    erfc(z.abs() / std::f64::consts::SQRT_2).clamp(0.0, 1.0)
}

pub(crate) fn chi_square_upper_tail(statistic: f64, degrees_of_freedom: usize) -> f64 {
    if !statistic.is_finite() || statistic < 0.0 || degrees_of_freedom == 0 {
        return f64::NAN;
    }
    let x = statistic / 2.0;
    if degrees_of_freedom.is_multiple_of(2) {
        let terms = degrees_of_freedom / 2;
        let mut term = 1.0;
        let mut sum = term;
        for index in 1..terms {
            term *= x / index as f64;
            sum += term;
        }
        (-x).exp() * sum
    } else {
        let mut result = erfc(x.sqrt());
        let mut shape = 0.5;
        let target = degrees_of_freedom as f64 / 2.0;
        let mut term = (-x).exp() * x.sqrt() / gamma_three_halves();
        while shape < target {
            result += term;
            shape += 1.0;
            term *= x / shape;
        }
        result.clamp(0.0, 1.0)
    }
}

fn gamma_three_halves() -> f64 {
    std::f64::consts::PI.sqrt() / 2.0
}

/// Complementary error function approximation with maximum error near 1e-7.
fn erfc(value: f64) -> f64 {
    if value.is_nan() {
        return f64::NAN;
    }
    if value < 0.0 {
        return 2.0 - erfc(-value);
    }
    let t = 1.0 / (1.0 + 0.5 * value);
    t * (-value * value - 1.265_512_23
        + t * (1.000_023_68
            + t * (0.374_091_96
                + t * (0.096_784_18
                    + t * (-0.186_288_06
                        + t * (0.278_868_07
                            + t * (-1.135_203_98
                                + t * (1.488_515_87 + t * (-0.822_152_23 + t * 0.170_872_77)))))))))
        .exp()
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

    #[test]
    fn average_ranks_handle_ties() {
        let (ranks, tie_sum) = average_ranks(&[3.0, 1.0, 1.0, 4.0]).unwrap();
        assert_eq!(ranks, vec![3.0, 1.5, 1.5, 4.0]);
        assert_eq!(tie_sum, 6.0);
    }

    #[test]
    fn distribution_tails_match_known_values() {
        assert!((standard_normal_two_sided_p(1.959_963_984_54) - 0.05).abs() < 2e-7);
        assert!((chi_square_upper_tail(3.841_458_820_69, 1) - 0.05).abs() < 2e-7);
        assert!((chi_square_upper_tail(5.991_464_547_11, 2) - 0.05).abs() < 2e-7);
        assert!((chi_square_upper_tail(7.814_727_903_25, 3) - 0.05).abs() < 2e-7);
    }
}
