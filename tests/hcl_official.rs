use pctsea::{AnyAtlas, AtlasData, HUMAN_CELL_LANDSCAPE};

/// Opt-in validation against the official download; ignored in normal CI.
///
/// Run with:
/// `PCTSEA_HCL_PATH=/path/HCL_Fig1_adata.h5ad cargo test --test hcl_official -- --ignored`
#[test]
#[ignore = "requires the 792 MiB official HCL download and paired cell-info workbook"]
fn official_hcl_opens_and_reads_expression() {
    let path = std::env::var_os("PCTSEA_HCL_PATH")
        .map(std::path::PathBuf::from)
        .expect("set PCTSEA_HCL_PATH to the official H5AD");
    let atlas = AnyAtlas::open(path).unwrap();
    let info = atlas.info();
    assert_eq!(info.cells, 599_926);
    assert_eq!(info.genes, 27_341);
    assert_eq!(info.cell_types, 63);
    assert_eq!(info.datasets, 106);
    assert!(
        atlas
            .gene_expression("CD3D")
            .unwrap()
            .is_some_and(|values| !values.is_empty())
    );
    assert_eq!(HUMAN_CELL_LANDSCAPE.expected_bytes, 830_846_460);
    assert_eq!(AtlasData::info(&atlas), info);
}
