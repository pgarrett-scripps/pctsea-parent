use crate::{AnalysisConfig, AnalysisResult, PlotConfig, enrichment_svg};
use std::fmt::Write;

/// Renders a complete, self-contained HTML analysis report.
pub fn html_report(
    result: &AnalysisResult,
    analysis_config: &AnalysisConfig,
    plot_config: &PlotConfig,
) -> String {
    let mut html = String::new();
    let plot = enrichment_svg(result, plot_config);
    let passing_rate = percent(result.cells_passing, result.cells_scored);
    let datasets = if analysis_config.datasets.is_empty() {
        "all datasets".to_string()
    } else {
        analysis_config.datasets.join(", ")
    };
    let threshold = analysis_config
        .score_threshold
        .map(format_number)
        .unwrap_or_else(|| "none".into());
    write!(
        html,
        r##"<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>pCTSEA analysis report</title>
<style>
:root {{ color-scheme: light; --ink: #172033; --muted: #526072; --line: #d8dee8; --soft: #f5f7fa; --positive: #b42318; --negative: #175cd3; --significant: #067647 }}
* {{ box-sizing: border-box }}
body {{ margin: 0; background: #eef2f6; color: var(--ink); font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; line-height: 1.45 }}
main {{ width: min(1180px, calc(100% - 32px)); margin: 32px auto 64px }}
header, section {{ background: #fff; border: 1px solid var(--line); border-radius: 12px; margin-bottom: 18px }}
header {{ padding: 28px 32px }}
section {{ padding: 24px }}
h1, h2 {{ margin: 0; letter-spacing: -0.02em }}
h1 {{ font-size: clamp(26px, 4vw, 38px) }}
h2 {{ font-size: 20px; margin-bottom: 16px }}
p {{ color: var(--muted); margin: 8px 0 0 }}
.summary {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; margin-top: 24px }}
.metric {{ background: var(--soft); border: 1px solid var(--line); border-radius: 8px; padding: 14px }}
.metric dt {{ color: var(--muted); font-size: 12px; margin-bottom: 3px }}
.metric dd {{ font-size: 22px; font-weight: 700; margin: 0 }}
.config {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap: 8px 20px; margin: 18px 0 0 }}
.config div {{ color: var(--muted); font-size: 13px }}
.config strong {{ color: var(--ink) }}
.plot {{ overflow-x: auto }}
.plot svg {{ display: block; height: auto; max-width: 100%; margin: 0 auto }}
.table-wrap {{ overflow: auto; max-height: 720px; border: 1px solid var(--line); border-radius: 8px }}
table {{ width: 100%; border-collapse: collapse; font-size: 13px; font-variant-numeric: tabular-nums }}
th {{ background: var(--soft); color: var(--muted); font-weight: 600; position: sticky; top: 0; text-align: left; z-index: 1 }}
th, td {{ border-bottom: 1px solid var(--line); padding: 9px 10px; white-space: nowrap }}
tbody tr:last-child td {{ border-bottom: 0 }}
tbody tr.significant {{ background: #ecfdf3 }}
.number {{ text-align: right }}
.positive {{ color: var(--positive); font-weight: 650 }}
.negative {{ color: var(--negative); font-weight: 650 }}
.q-significant {{ color: var(--significant); font-weight: 750 }}
details {{ margin-top: 16px }}
summary {{ color: var(--ink); cursor: pointer; font-weight: 650; margin-bottom: 12px }}
.genes {{ white-space: normal; min-width: 240px }}
.note {{ border-left: 3px solid #f79009; padding-left: 12px }}
footer {{ color: var(--muted); font-size: 12px; padding: 8px 4px }}
@media print {{ body {{ background: #fff }} main {{ margin: 0; width: 100% }} header, section {{ break-inside: avoid; border-color: #bbb }} .table-wrap {{ max-height: none; overflow: visible }} th {{ position: static }} }}
</style>
</head>
<body>
<main>
<header>
<h1>pCTSEA analysis report</h1>
<p>Proteomic cell-type enrichment with ranked-cell and score-distribution statistics.</p>
<dl class="summary">
<div class="metric"><dt>Matched genes</dt><dd>{}</dd></div>
<div class="metric"><dt>Missing genes</dt><dd>{}</dd></div>
<div class="metric"><dt>Cells scored</dt><dd>{}</dd></div>
<div class="metric"><dt>Cells passing</dt><dd>{} ({passing_rate})</dd></div>
<div class="metric"><dt>Cell types</dt><dd>{}</dd></div>
</dl>
<div class="config">
<div><strong>Scoring:</strong> {}</div>
<div><strong>Minimum genes per cell:</strong> {}</div>
<div><strong>Score threshold:</strong> {threshold}</div>
<div><strong>Permutations:</strong> {}</div>
<div><strong>Random seed:</strong> {}</div>
<div><strong>Datasets:</strong> {}</div>
<div><strong>Positive pairs only:</strong> {}</div>
</div>
</header>
<section>
<h2>Enrichment overview</h2>
<div class="plot">{plot}</div>
</section>
"##,
        result.matched_genes.len(),
        result.missing_genes.len(),
        result.cells_scored,
        result.cells_passing,
        result.cell_types.len(),
        analysis_config.method,
        analysis_config.min_genes_per_cell,
        analysis_config.permutations,
        analysis_config.random_seed,
        escape_html(&datasets),
        if analysis_config.positive_pairs_only {
            "yes"
        } else {
            "no"
        }
    )
    .unwrap();

    html.push_str(
        r##"<section>
<h2>Cell-type results</h2>
<p>Positive normalized enrichment scores indicate concentration near the high-scoring end of the ranking. Rows with permutation FDR at or below the configured plot threshold are highlighted.</p>
<div class="table-wrap"><table>
<thead><tr><th>Cell type</th><th class="number">Total</th><th class="number">Passing</th><th class="number">Pass rate</th><th class="number">NES</th><th class="number">Permutation p</th><th class="number">Permutation FDR</th><th class="number">Hypergeometric FDR</th><th>Top genes</th></tr></thead>
<tbody>
"##,
    );
    for row in &result.cell_types {
        let significant =
            row.permutation_fdr.is_finite() && row.permutation_fdr <= plot_config.fdr_threshold;
        let direction = if row.normalized_enrichment_score >= 0.0 {
            "positive"
        } else {
            "negative"
        };
        let top_genes = row
            .top_genes
            .iter()
            .map(|(gene, count)| format!("{} ({count})", escape_html(gene)))
            .collect::<Vec<_>>()
            .join(", ");
        writeln!(
            html,
            r##"<tr{}><td>{}</td><td class="number">{}</td><td class="number">{}</td><td class="number">{}</td><td class="number {direction}">{}</td><td class="number">{}</td><td class="number{}">{}</td><td class="number">{}</td><td class="genes">{}</td></tr>"##,
            if significant {
                " class=\"significant\""
            } else {
                ""
            },
            escape_html(&row.cell_type),
            row.cells_total,
            row.cells_passing,
            percent(row.cells_passing, row.cells_total),
            format_number(row.normalized_enrichment_score),
            format_probability(row.permutation_p_value),
            if significant { " q-significant" } else { "" },
            format_probability(row.permutation_fdr),
            format_probability(row.hypergeometric_fdr),
            top_genes
        )
        .unwrap();
    }
    html.push_str("</tbody></table></div>\n</section>\n");

    write!(
        html,
        r##"<section>
<h2>Score-distribution follow-up</h2>
<p>Kruskal-Wallis H = {}, df = {}, p = {}. Dunn comparisons use combined ranks and Benjamini-Hochberg correction across all {} pairs.</p>
<details open><summary>Pairwise comparisons</summary>
<div class="table-wrap"><table>
<thead><tr><th>Cell type A</th><th>Cell type B</th><th class="number">n A</th><th class="number">n B</th><th class="number">Median A</th><th class="number">Median B</th><th class="number">Median difference</th><th class="number">Dunn z</th><th class="number">p</th><th class="number">FDR</th></tr></thead>
<tbody>
"##,
        format_number(result.score_distribution_test.statistic),
        result.score_distribution_test.degrees_of_freedom,
        format_probability(result.score_distribution_test.p_value),
        result.post_hoc_comparisons.len()
    )
    .unwrap();
    for comparison in &result.post_hoc_comparisons {
        let significant = comparison.fdr.is_finite() && comparison.fdr <= plot_config.fdr_threshold;
        let direction = if comparison.median_difference >= 0.0 {
            "positive"
        } else {
            "negative"
        };
        writeln!(
            html,
            r##"<tr{}><td>{}</td><td>{}</td><td class="number">{}</td><td class="number">{}</td><td class="number">{}</td><td class="number">{}</td><td class="number {direction}">{}</td><td class="number">{}</td><td class="number">{}</td><td class="number{}">{}</td></tr>"##,
            if significant {
                " class=\"significant\""
            } else {
                ""
            },
            escape_html(&comparison.cell_type_a),
            escape_html(&comparison.cell_type_b),
            comparison.cells_a,
            comparison.cells_b,
            format_number(comparison.median_score_a),
            format_number(comparison.median_score_b),
            format_number(comparison.median_difference),
            format_number(comparison.z_score),
            format_probability(comparison.p_value),
            if significant { " q-significant" } else { "" },
            format_probability(comparison.fdr)
        )
        .unwrap();
    }
    html.push_str("</tbody></table></div>\n</details>\n</section>\n");

    write!(
        html,
        r##"<section>
<h2>Query coverage</h2>
<details><summary>Matched genes ({})</summary><p>{}</p></details>
<details><summary>Missing genes ({})</summary><p>{}</p></details>
</section>
<section>
<h2>Interpretation notes</h2>
<p>The enrichment tests and the post-hoc rank tests answer different questions. Use the enrichment table for the primary cell-type result. Use Dunn comparisons to understand which score distributions differ after the omnibus test.</p>
<p class="note">The score-distribution tests treat cells as independent observations. Donor-level or batch-level studies need replicate-aware validation before publication.</p>
</section>
<footer>Generated by pCTSEA. This report is self-contained and can be archived or shared as one HTML file.</footer>
</main>
</body>
</html>
"##,
        result.matched_genes.len(),
        escaped_list(&result.matched_genes),
        result.missing_genes.len(),
        escaped_list(&result.missing_genes)
    )
    .unwrap();
    html
}

fn escaped_list(values: &[String]) -> String {
    if values.is_empty() {
        "None".into()
    } else {
        values
            .iter()
            .map(|value| escape_html(value))
            .collect::<Vec<_>>()
            .join(", ")
    }
}

fn percent(numerator: usize, denominator: usize) -> String {
    if denominator == 0 {
        "NA".into()
    } else {
        format!("{:.1}%", numerator as f64 * 100.0 / denominator as f64)
    }
}

fn format_number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.4}")
    } else {
        "NA".into()
    }
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

fn escape_html(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escapes_html_text() {
        assert_eq!(
            escape_html("T & <NK> \"cells\""),
            "T &amp; &lt;NK&gt; &quot;cells&quot;"
        );
    }

    #[test]
    fn percentages_handle_empty_groups() {
        assert_eq!(percent(1, 4), "25.0%");
        assert_eq!(percent(0, 0), "NA");
    }
}
