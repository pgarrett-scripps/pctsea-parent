use crate::AnalysisResult;
use std::fmt::Write;

/// Controls the standalone SVG enrichment plot.
#[derive(Clone, Debug)]
pub struct PlotConfig {
    pub width: usize,
    /// Zero includes every cell type.
    pub max_cell_types: usize,
    pub fdr_threshold: f64,
}

impl Default for PlotConfig {
    fn default() -> Self {
        Self {
            width: 1_000,
            max_cell_types: 30,
            fdr_threshold: 0.05,
        }
    }
}

/// Renders a dependency-free SVG summary of normalized enrichment scores.
pub fn enrichment_svg(result: &AnalysisResult, config: &PlotConfig) -> String {
    let width = config.width.max(640);
    let limit = if config.max_cell_types == 0 {
        usize::MAX
    } else {
        config.max_cell_types
    };
    let mut rows: Vec<_> = result
        .cell_types
        .iter()
        .filter(|row| row.normalized_enrichment_score.is_finite())
        .take(limit)
        .collect();
    rows.sort_by(|left, right| {
        left.normalized_enrichment_score
            .total_cmp(&right.normalized_enrichment_score)
    });

    let top = 112.0;
    let row_height = 30.0;
    let bottom = 64.0;
    let height = (top + rows.len().max(1) as f64 * row_height + bottom).ceil() as usize;
    let left = (width as f64 * 0.28).clamp(190.0, 300.0);
    let right = 130.0;
    let plot_width = width as f64 - left - right;
    let maximum = rows
        .iter()
        .map(|row| row.normalized_enrichment_score.abs())
        .fold(0.0_f64, f64::max)
        .max(0.1);
    let scale = |value: f64| left + (value / maximum + 1.0) * plot_width / 2.0;
    let zero = scale(0.0);
    let mut svg = String::new();
    writeln!(
        svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title description">"##
    )
    .unwrap();
    writeln!(
        svg,
        r##"<title id="title">pCTSEA cell-type enrichment</title>"##
    )
    .unwrap();
    writeln!(
        svg,
        r##"<desc id="description">Normalized enrichment scores by cell type. Saturated marks have permutation FDR at or below the configured threshold.</desc>"##
    )
    .unwrap();
    writeln!(
        svg,
        r##"<rect width="100%" height="100%" fill="#ffffff"/>"##
    )
    .unwrap();
    writeln!(
        svg,
        r##"<text x="28" y="36" fill="#172033" font-family="system-ui, sans-serif" font-size="22" font-weight="700">pCTSEA cell-type enrichment</text>"##
    )
    .unwrap();
    writeln!(
        svg,
        r##"<text x="28" y="62" fill="#526072" font-family="system-ui, sans-serif" font-size="13">Normalized enrichment score with permutation FDR. Kruskal-Wallis p = {}</text>"##,
        format_probability(result.score_distribution_test.p_value)
    )
    .unwrap();
    writeln!(
        svg,
        r##"<text x="{}" y="91" text-anchor="end" fill="#526072" font-family="system-ui, sans-serif" font-size="12">cell type</text>"##,
        left - 14.0
    )
    .unwrap();
    writeln!(
        svg,
        r##"<text x="{}" y="91" fill="#526072" font-family="system-ui, sans-serif" font-size="12">permutation q</text>"##,
        width as f64 - right + 18.0
    )
    .unwrap();

    for value in [-maximum, 0.0, maximum] {
        let x = scale(value);
        writeln!(
            svg,
            r##"<line x1="{x:.2}" y1="100" x2="{x:.2}" y2="{}" stroke="{}" stroke-width="1"/>"##,
            height as f64 - bottom + 5.0,
            if value == 0.0 { "#718096" } else { "#d8dee8" }
        )
        .unwrap();
        writeln!(
            svg,
            r##"<text x="{x:.2}" y="{}" text-anchor="middle" fill="#526072" font-family="system-ui, sans-serif" font-size="12">{value:.2}</text>"##,
            height as f64 - 31.0
        )
        .unwrap();
    }

    if rows.is_empty() {
        writeln!(
            svg,
            r##"<text x="{}" y="{}" text-anchor="middle" fill="#526072" font-family="system-ui, sans-serif" font-size="14">No finite enrichment scores to plot</text>"##,
            width / 2,
            height / 2
        )
        .unwrap();
    }
    for (index, row) in rows.iter().enumerate() {
        let y = top + index as f64 * row_height + row_height / 2.0;
        let x = scale(row.normalized_enrichment_score);
        let significant =
            row.permutation_fdr.is_finite() && row.permutation_fdr <= config.fdr_threshold;
        let color = if row.normalized_enrichment_score >= 0.0 {
            if significant { "#b42318" } else { "#d7a29d" }
        } else if significant {
            "#175cd3"
        } else {
            "#9bb7df"
        };
        writeln!(
            svg,
            r##"<text x="{}" y="{:.2}" text-anchor="end" dominant-baseline="middle" fill="#172033" font-family="system-ui, sans-serif" font-size="12">{}</text>"##,
            left - 14.0,
            y,
            escape_xml(&row.cell_type)
        )
        .unwrap();
        writeln!(
            svg,
            r##"<line x1="{zero:.2}" y1="{y:.2}" x2="{x:.2}" y2="{y:.2}" stroke="{color}" stroke-width="5" stroke-linecap="round"/>"##
        )
        .unwrap();
        writeln!(
            svg,
            r##"<circle cx="{x:.2}" cy="{y:.2}" r="5" fill="{color}"/>"##
        )
        .unwrap();
        writeln!(
            svg,
            r##"<text x="{}" y="{y:.2}" dominant-baseline="middle" fill="#172033" font-family="ui-monospace, monospace" font-size="12">{}</text>"##,
            width as f64 - right + 18.0,
            format_probability(row.permutation_fdr)
        )
        .unwrap();
    }
    writeln!(
        svg,
        r##"<text x="{}" y="{}" text-anchor="middle" fill="#172033" font-family="system-ui, sans-serif" font-size="13">normalized enrichment score</text>"##,
        left + plot_width / 2.0,
        height - 8
    )
    .unwrap();
    svg.push_str("</svg>\n");
    svg
}

fn format_probability(value: f64) -> String {
    if !value.is_finite() {
        "NA".into()
    } else if value < 0.0001 {
        format!("{value:.2e}")
    } else {
        format!("{value:.4}")
    }
}

fn escape_xml(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CellTypeResult, KruskalWallisResult};

    #[test]
    fn svg_contains_labels_and_escapes_xml() {
        let result = AnalysisResult {
            matched_genes: Vec::new(),
            missing_genes: Vec::new(),
            cells_considered: 2,
            cells_scored: 2,
            cells_passing: 2,
            cell_scores: Vec::new(),
            cell_types: vec![CellTypeResult {
                cell_type: "T & NK".into(),
                cells_total: 2,
                cells_passing: 2,
                hypergeometric_p_value: 1.0,
                hypergeometric_fdr: 1.0,
                enrichment_score: 0.5,
                normalized_enrichment_score: 1.25,
                permutation_p_value: 0.01,
                permutation_fdr: 0.02,
                top_genes: Vec::new(),
            }],
            score_distribution_test: KruskalWallisResult {
                statistic: 4.0,
                degrees_of_freedom: 1,
                p_value: 0.0455,
            },
            post_hoc_comparisons: Vec::new(),
        };
        let svg = enrichment_svg(&result, &PlotConfig::default());
        assert!(svg.starts_with("<svg"));
        assert!(svg.contains("T &amp; NK"));
        assert!(svg.contains("1.25"));
        assert!(svg.ends_with("</svg>\n"));
    }
}
