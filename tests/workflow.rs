use pctsea::{AnalysisConfig, Atlas, AtlasBuilder, GeneQuery, analyze};
use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

fn temp_file(name: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("pctsea-workflow-{nonce}-{name}"))
}

#[test]
fn public_api_runs_end_to_end() {
    let cells = temp_file("cells.tsv");
    let expressions = temp_file("expressions.tsv");
    let atlas_path = temp_file("atlas.pctsea");
    fs::write(
        &cells,
        "cell\tcell_type\tdataset\nt1\tT\tdemo\nt2\tT\tdemo\nb1\tB\tdemo\nb2\tB\tdemo\n",
    )
    .unwrap();
    fs::write(
        &expressions,
        concat!(
            "gene\tcell\texpression\n",
            "G1\tt1\t4\nG2\tt1\t3\nG3\tt1\t2\n",
            "G1\tt2\t3\nG2\tt2\t2\nG3\tt2\t1\n",
            "G1\tb1\t1\nG2\tb1\t2\nG3\tb1\t4\n",
            "G1\tb2\t1\nG2\tb2\t3\nG3\tb2\t4\n",
        ),
    )
    .unwrap();
    AtlasBuilder::pack_tsv(&cells, &expressions, &atlas_path).unwrap();

    let atlas = Atlas::open(&atlas_path).unwrap();
    let query = GeneQuery::new([("G1", 4.0), ("G2", 3.0), ("G3", 1.0)]).unwrap();
    let config = AnalysisConfig {
        min_genes_per_cell: 3,
        score_threshold: None,
        permutations: 100,
        random_seed: 7,
        ..AnalysisConfig::default()
    };
    let result = analyze(&atlas, &query, &config).unwrap();

    assert_eq!(result.cells_scored, 4);
    assert_eq!(result.cells_passing, 4);
    let t_cells = result
        .cell_types
        .iter()
        .find(|result| result.cell_type == "T")
        .unwrap();
    let b_cells = result
        .cell_types
        .iter()
        .find(|result| result.cell_type == "B")
        .unwrap();
    assert!(t_cells.enrichment_score > 0.0);
    assert!(b_cells.enrichment_score < 0.0);
    assert!(t_cells.permutation_p_value.is_finite());

    for path in [cells, expressions, atlas_path] {
        fs::remove_file(path).unwrap();
    }
}
