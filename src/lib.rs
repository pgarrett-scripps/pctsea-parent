//! A self-contained, file-backed cell-type enrichment library.
//!
//! The native atlas is sparse and indexed by gene. Opening an [`Atlas`] reads
//! metadata and the gene directory only; [`Atlas::load_query`] seeks directly
//! to expression postings for the requested genes.

mod analysis;
mod atlas;
mod catalog;
mod error;
mod h5ad;
mod input;
mod stats;

pub use analysis::{
    AnalysisConfig, AnalysisResult, CellScore, CellTypeResult, ScoringMethod, analyze,
};
pub use atlas::{Atlas, AtlasBuilder, AtlasData, AtlasInfo, Cell, GeneExpression, LoadedQuery};
pub use catalog::{
    AtlasDownload, AtlasFileKind, AtlasRelease, HUMAN_CELL_LANDSCAPE,
    HUMAN_CELL_LANDSCAPE_ANNOTATIONS, default_atlas_path, download_atlas,
};
pub use error::{Error, Result};
pub use h5ad::{AnyAtlas, H5adAtlas};
pub use input::{GeneQuery, QueryGene};
