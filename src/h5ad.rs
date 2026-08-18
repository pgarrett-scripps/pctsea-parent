use crate::{
    Atlas, AtlasData, AtlasInfo, Cell, Error, GeneExpression, GeneQuery, LoadedQuery, Result,
};
use anndata::data::{DynArray, DynCscMatrix, DynCsrMatrix, SelectInfoElem};
use anndata::{AnnData, AnnDataOp, ArrayData, ArrayElemOp, Backend};
use anndata_hdf5::H5;
use calamine::{DataType as CalamineDataType, Reader, open_workbook_auto};
use hdf5::types::{FixedAscii, VarLenUnicode};
use hdf5::{Dataset, File as Hdf5File};
use polars::prelude::{DataFrame, DataType};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

const NATIVE_MAGIC: &[u8; 8] = b"PCTSEA\0\x01";
const HDF5_MAGIC: &[u8; 8] = b"\x89HDF\r\n\x1a\n";
const HCL_ANNOTATIONS_FILE: &str = "HCL_Fig1_cell_Info.xlsx";
const METADATA_MAGIC: &[u8; 8] = b"PCTMETA1";
const MAX_CELL_TYPE_BYTES: usize = 1024 * 1024;

#[derive(Clone, hdf5::H5Type)]
#[repr(C, packed)]
struct LegacyObs {
    index: FixedAscii<44>,
    batch: u8,
    tissue: u8,
    n_genes: u64,
    n_counts: f32,
}

#[derive(Clone, hdf5::H5Type)]
#[repr(C, packed)]
struct LegacyVar {
    index: FixedAscii<22>,
    n_cells: u64,
}

/// An atlas opened from either H5AD or the compact native cache format.
#[derive(Debug)]
pub enum AnyAtlas {
    H5ad(H5adAtlas),
    Native(Atlas),
}

impl AnyAtlas {
    /// Detects the format from magic bytes rather than relying on its extension.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let mut file = File::open(path)?;
        let mut magic = [0_u8; 8];
        file.read_exact(&mut magic)?;
        match &magic {
            NATIVE_MAGIC => Atlas::open(path).map(Self::Native),
            HDF5_MAGIC => H5adAtlas::open(path).map(Self::H5ad),
            _ => Err(Error::InvalidAtlas(format!(
                "{} is neither a PCTSEA atlas nor an HDF5/H5AD file",
                path.display()
            ))),
        }
    }

    pub fn cells(&self) -> &[Cell] {
        AtlasData::cells(self)
    }

    pub fn info(&self) -> AtlasInfo {
        AtlasData::info(self)
    }

    pub fn gene_expression(&self, gene: &str) -> Result<Option<Vec<GeneExpression>>> {
        AtlasData::gene_expression(self, gene)
    }
}

impl AtlasData for AnyAtlas {
    fn cells(&self) -> &[Cell] {
        match self {
            Self::H5ad(atlas) => atlas.cells(),
            Self::Native(atlas) => atlas.cells(),
        }
    }

    fn info(&self) -> AtlasInfo {
        match self {
            Self::H5ad(atlas) => atlas.info(),
            Self::Native(atlas) => atlas.info(),
        }
    }

    fn gene_expression(&self, gene: &str) -> Result<Option<Vec<GeneExpression>>> {
        match self {
            Self::H5ad(atlas) => atlas.gene_expression(gene),
            Self::Native(atlas) => atlas.gene_expression(gene),
        }
    }

    fn load_query(&self, query: &GeneQuery) -> Result<LoadedQuery> {
        match self {
            Self::H5ad(atlas) => atlas.load_query(query),
            Self::Native(atlas) => atlas.load_query(query),
        }
    }
}

/// A read-only, file-backed H5AD atlas.
///
/// Cell annotations and gene names are resident. Expression values remain in
/// HDF5 and only columns needed by a query are read.
#[derive(Debug)]
pub struct H5adAtlas {
    storage: H5adStorage,
    cells: Vec<Cell>,
    genes: Vec<String>,
    gene_lookup: HashMap<String, usize>,
}

#[derive(Debug)]
enum H5adStorage {
    Modern(AnnData<H5>),
    LegacyDense(Dataset),
}

impl H5adAtlas {
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        match Self::open_modern(path) {
            Ok(atlas) => Ok(atlas),
            Err(modern_error) => {
                let annotations = path.with_file_name(HCL_ANNOTATIONS_FILE);
                Self::open_legacy(path, &annotations).map_err(|legacy_error| {
                    Error::InvalidAtlas(format!(
                        "could not read modern H5AD ({modern_error}); legacy H5AD fallback also failed ({legacy_error})"
                    ))
                })
            }
        }
    }

    /// Opens a legacy H5AD with an explicit cell-annotation workbook sidecar.
    pub fn open_with_annotations(
        path: impl AsRef<Path>,
        annotations: impl AsRef<Path>,
    ) -> Result<Self> {
        Self::open_modern(path.as_ref())
            .or_else(|_| Self::open_legacy(path.as_ref(), annotations.as_ref()))
    }

    fn open_modern(path: &Path) -> Result<Self> {
        let store = H5::open(path).map_err(h5ad_error)?;
        let data = AnnData::<H5>::open(store).map_err(h5ad_error)?;
        if data.x().is_none() {
            return Err(Error::InvalidAtlas("H5AD file has no X matrix".into()));
        }
        let names = data.obs_names().into_vec();
        let obs = data.read_obs().map_err(h5ad_error)?;
        let cell_types = required_string_column(
            &obs,
            &[
                "cell_type",
                "celltype",
                "cell_types",
                "cell_ontology_class",
                "annotation",
                "cluster",
            ],
            "cell type",
        )?;
        let datasets = optional_string_column(
            &obs,
            &["dataset", "tissue", "organ", "sample", "batch", "study"],
        )
        .unwrap_or_else(|| vec![String::new(); names.len()]);
        if names.len() != cell_types.len() || names.len() != datasets.len() {
            return Err(Error::InvalidAtlas(
                "H5AD observation annotations do not match the X matrix".into(),
            ));
        }
        let cells = names
            .into_iter()
            .zip(cell_types)
            .zip(datasets)
            .map(|((name, cell_type), dataset)| Cell {
                name,
                cell_type,
                dataset,
            })
            .collect();
        let genes = data.var_names().into_vec();
        let mut gene_lookup = HashMap::with_capacity(genes.len());
        for (index, gene) in genes.iter().enumerate() {
            gene_lookup
                .entry(gene.trim().to_uppercase())
                .or_insert(index);
        }
        Ok(Self {
            storage: H5adStorage::Modern(data),
            cells,
            genes,
            gene_lookup,
        })
    }

    fn open_legacy(path: &Path, annotations: &Path) -> Result<Self> {
        if !annotations.is_file() {
            return Err(Error::InvalidAtlas(format!(
                "legacy HCL H5AD requires annotation sidecar {}; run `pctsea atlas download hcl` or use H5adAtlas::open_with_annotations",
                annotations.display()
            )));
        }
        let file = Hdf5File::open(path).map_err(h5ad_error)?;
        let x = file.dataset("X").map_err(h5ad_error)?;
        let shape = x.shape();
        if shape.len() != 2 {
            return Err(Error::InvalidAtlas("legacy H5AD X is not a matrix".into()));
        }
        let observations = file
            .dataset("obs")
            .map_err(h5ad_error)?
            .read_1d::<LegacyObs>()
            .map_err(h5ad_error)?;
        let variables = file
            .dataset("var")
            .map_err(h5ad_error)?
            .read_1d::<LegacyVar>()
            .map_err(h5ad_error)?;
        if shape != [observations.len(), variables.len()] {
            return Err(Error::InvalidAtlas(
                "legacy H5AD obs/var dimensions do not match X".into(),
            ));
        }
        let batch_categories = read_string_dataset(&file, "uns/batch_categories")?;
        let tissue_categories = read_string_dataset(&file, "uns/tissue_categories")?;
        let cell_types = match read_metadata_cache(path, annotations, observations.len()) {
            Ok(Some(cell_types)) => cell_types,
            _ => {
                let mut annotation_map = read_hcl_annotations(annotations)?;
                let mut missing = Vec::new();
                let cell_types: Vec<String> = observations
                    .iter()
                    .map(|observation| {
                        let name = observation.index.as_str();
                        annotation_map
                            .remove(legacy_cell_key(name))
                            .unwrap_or_else(|| {
                                if missing.len() < 5 {
                                    missing.push(name.to_owned());
                                }
                                "Unknown".into()
                            })
                    })
                    .collect();
                if !missing.is_empty() {
                    return Err(Error::InvalidAtlas(format!(
                        "annotation sidecar does not contain H5AD cells such as {}",
                        missing.join(", ")
                    )));
                }
                let _ = write_metadata_cache(path, annotations, &cell_types);
                cell_types
            }
        };
        let mut cells = Vec::with_capacity(observations.len());
        for (observation, cell_type) in observations.into_iter().zip(cell_types) {
            let name = observation.index.as_str().to_owned();
            let batch = batch_categories
                .get(observation.batch as usize)
                .ok_or_else(|| Error::InvalidAtlas("legacy H5AD batch code is invalid".into()))?;
            let tissue = tissue_categories
                .get(observation.tissue as usize)
                .ok_or_else(|| Error::InvalidAtlas("legacy H5AD tissue code is invalid".into()))?;
            cells.push(Cell {
                name,
                cell_type,
                dataset: format!("{tissue}/{batch}"),
            });
        }
        let genes: Vec<String> = variables
            .iter()
            .map(|variable| variable.index.as_str().to_owned())
            .collect();
        let gene_lookup = gene_lookup(&genes);
        Ok(Self {
            storage: H5adStorage::LegacyDense(x),
            cells,
            genes,
            gene_lookup,
        })
    }

    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    pub fn info(&self) -> AtlasInfo {
        AtlasInfo {
            format: "h5ad".into(),
            format_version: None,
            cells: self.cells.len(),
            genes: self.genes.len(),
            non_zero_expressions: None,
            cell_types: self
                .cells
                .iter()
                .map(|cell| &cell.cell_type)
                .collect::<HashSet<_>>()
                .len(),
            datasets: self
                .cells
                .iter()
                .map(|cell| &cell.dataset)
                .filter(|dataset| !dataset.is_empty())
                .collect::<HashSet<_>>()
                .len(),
        }
    }

    pub fn gene_expression(&self, gene: &str) -> Result<Option<Vec<GeneExpression>>> {
        if !self.gene_lookup.contains_key(&gene.trim().to_uppercase()) {
            return Ok(None);
        }
        let query = GeneQuery::new([(gene, 1.0)])?;
        let loaded = self.load_query(&query)?;
        Ok(Some(
            loaded
                .cell_values
                .into_iter()
                .enumerate()
                .filter_map(|(cell_index, row)| {
                    let value = row[0];
                    (value != 0.0).then_some(GeneExpression { cell_index, value })
                })
                .collect(),
        ))
    }

    pub fn load_query(&self, query: &GeneQuery) -> Result<LoadedQuery> {
        let mut gene_names = Vec::new();
        let mut query_values = Vec::new();
        let mut indices = Vec::new();
        let mut missing_genes = Vec::new();
        for gene in &query.genes {
            if let Some(&index) = self.gene_lookup.get(&gene.name) {
                gene_names.push(gene.name.clone());
                query_values.push(gene.value);
                indices.push(index);
            } else {
                missing_genes.push(gene.name.clone());
            }
        }
        if indices.is_empty() {
            return Err(Error::InvalidConfig(
                "none of the query genes occur in the atlas".into(),
            ));
        }
        let cell_values = match &self.storage {
            H5adStorage::Modern(data) => {
                let selection = [SelectInfoElem::from(..), SelectInfoElem::from(indices)];
                let matrix = data
                    .x()
                    .slice::<ArrayData, _>(&selection)
                    .map_err(h5ad_error)?
                    .ok_or_else(|| Error::InvalidAtlas("H5AD file has no X matrix".into()))?;
                matrix_to_rows(matrix, self.cells.len(), gene_names.len())?
            }
            H5adStorage::LegacyDense(dataset) => {
                read_legacy_columns(dataset, &indices, self.cells.len())?
            }
        };
        Ok(LoadedQuery {
            gene_names,
            query_values,
            cell_values,
            missing_genes,
        })
    }
}

impl AtlasData for H5adAtlas {
    fn cells(&self) -> &[Cell] {
        self.cells()
    }

    fn info(&self) -> AtlasInfo {
        self.info()
    }

    fn gene_expression(&self, gene: &str) -> Result<Option<Vec<GeneExpression>>> {
        self.gene_expression(gene)
    }

    fn load_query(&self, query: &GeneQuery) -> Result<LoadedQuery> {
        self.load_query(query)
    }
}

fn gene_lookup(genes: &[String]) -> HashMap<String, usize> {
    let mut lookup = HashMap::with_capacity(genes.len());
    for (index, gene) in genes.iter().enumerate() {
        lookup.entry(gene.trim().to_uppercase()).or_insert(index);
    }
    lookup
}

fn metadata_cache_path(annotations: &Path) -> PathBuf {
    let mut path = annotations.as_os_str().to_os_string();
    path.push(".pctsea-meta");
    PathBuf::from(path)
}

fn read_metadata_cache(
    atlas: &Path,
    annotations: &Path,
    expected_cells: usize,
) -> Result<Option<Vec<String>>> {
    let cache = metadata_cache_path(annotations);
    if !cache.is_file() {
        return Ok(None);
    }
    let mut reader = BufReader::new(File::open(cache)?);
    let mut magic = [0_u8; 8];
    reader.read_exact(&mut magic)?;
    if &magic != METADATA_MAGIC {
        return Ok(None);
    }
    let atlas_bytes = read_cache_u64(&mut reader)?;
    let annotation_bytes = read_cache_u64(&mut reader)?;
    let cell_count = read_cache_u64(&mut reader)? as usize;
    if atlas_bytes != atlas.metadata()?.len()
        || annotation_bytes != annotations.metadata()?.len()
        || cell_count != expected_cells
    {
        return Ok(None);
    }
    let mut cell_types = Vec::with_capacity(cell_count);
    for _ in 0..cell_count {
        let length = read_cache_u32(&mut reader)? as usize;
        if length > MAX_CELL_TYPE_BYTES {
            return Ok(None);
        }
        let mut bytes = vec![0_u8; length];
        reader.read_exact(&mut bytes)?;
        let Ok(value) = String::from_utf8(bytes) else {
            return Ok(None);
        };
        cell_types.push(value);
    }
    Ok(Some(cell_types))
}

fn write_metadata_cache(atlas: &Path, annotations: &Path, cell_types: &[String]) -> Result<()> {
    let cache = metadata_cache_path(annotations);
    let mut part_name = cache.as_os_str().to_os_string();
    part_name.push(".part");
    let part = PathBuf::from(part_name);
    let mut writer = BufWriter::new(File::create(&part)?);
    writer.write_all(METADATA_MAGIC)?;
    writer.write_all(&atlas.metadata()?.len().to_le_bytes())?;
    writer.write_all(&annotations.metadata()?.len().to_le_bytes())?;
    writer.write_all(&(cell_types.len() as u64).to_le_bytes())?;
    for cell_type in cell_types {
        if cell_type.len() > u32::MAX as usize {
            return Err(Error::InvalidAtlas("cell type is too long to cache".into()));
        }
        writer.write_all(&(cell_type.len() as u32).to_le_bytes())?;
        writer.write_all(cell_type.as_bytes())?;
    }
    writer.flush()?;
    drop(writer);
    if cache.exists() {
        fs::remove_file(&cache)?;
    }
    fs::rename(part, cache)?;
    Ok(())
}

fn read_cache_u32(reader: &mut impl Read) -> Result<u32> {
    let mut bytes = [0_u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_cache_u64(reader: &mut impl Read) -> Result<u64> {
    let mut bytes = [0_u8; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

fn read_string_dataset(file: &Hdf5File, path: &str) -> Result<Vec<String>> {
    file.dataset(path)
        .map_err(h5ad_error)?
        .read_1d::<VarLenUnicode>()
        .map_err(h5ad_error)
        .map(|values| {
            values
                .iter()
                .map(|value| value.as_str().to_owned())
                .collect()
        })
}

fn read_hcl_annotations(path: &Path) -> Result<HashMap<String, String>> {
    let mut workbook = open_workbook_auto(path).map_err(h5ad_error)?;
    let mut annotations = HashMap::new();
    for sheet in workbook.sheet_names().to_vec() {
        let range = workbook.worksheet_range(&sheet).map_err(h5ad_error)?;
        let mut rows = range.rows();
        let headers = rows
            .next()
            .ok_or_else(|| Error::InvalidAtlas("annotation workbook sheet is empty".into()))?;
        let cell_index = headers
            .iter()
            .position(|name| {
                name.as_string().is_some_and(|name| {
                    matches!(
                        normalize_column(&name).as_str(),
                        "cellid" | "cell" | "cellnames" | "cellbarcode" | "barcode"
                    )
                })
            })
            .ok_or_else(|| {
                Error::InvalidAtlas("annotation workbook has no cell ID column".into())
            })?;
        let type_index = headers
            .iter()
            .position(|name| {
                name.as_string()
                    .is_some_and(|name| normalize_column(&name) == "celltype")
            })
            .ok_or_else(|| {
                Error::InvalidAtlas("annotation workbook has no cell type column".into())
            })?;
        for row in rows {
            let cell = row
                .get(cell_index)
                .and_then(CalamineDataType::as_string)
                .unwrap_or_default();
            let cell_type = row
                .get(type_index)
                .and_then(CalamineDataType::as_string)
                .unwrap_or_default();
            let cell = cell.trim();
            let cell_type = cell_type.trim();
            if cell.is_empty() || cell_type.is_empty() {
                return Err(Error::InvalidAtlas(
                    "annotation workbook contains an empty cell ID or cell type".into(),
                ));
            }
            annotations
                .entry(cell.to_owned())
                .or_insert_with(|| cell_type.to_owned());
        }
    }
    if annotations.is_empty() {
        return Err(Error::InvalidAtlas(
            "annotation workbook contains no annotation records".into(),
        ));
    }
    Ok(annotations)
}

fn legacy_cell_key(name: &str) -> &str {
    name.rsplit_once('-').map_or(name, |(prefix, suffix)| {
        if suffix.chars().all(|character| character.is_ascii_digit()) {
            prefix
        } else {
            name
        }
    })
}

fn read_legacy_columns(dataset: &Dataset, indices: &[usize], rows: usize) -> Result<Vec<Vec<f32>>> {
    let shape = dataset.shape();
    if shape.len() != 2 || shape[0] != rows {
        return Err(Error::InvalidAtlas(
            "legacy H5AD X has an unexpected shape".into(),
        ));
    }
    let columns = shape[1];
    if indices.iter().any(|&index| index >= columns) {
        return Err(Error::InvalidAtlas(
            "legacy H5AD query refers to a missing gene column".into(),
        ));
    }
    let mut output = vec![vec![0.0; indices.len()]; rows];
    let chunk_columns = dataset
        .chunk()
        .and_then(|shape| shape.get(1).copied())
        .unwrap_or(1)
        .max(1);
    let groups = legacy_column_groups(indices, chunk_columns);
    for (chunk_index, selected_columns) in groups {
        let source_start = chunk_index * chunk_columns;
        let source_end = (source_start + chunk_columns).min(columns);
        let values = dataset
            .read_slice_2d::<f32, _>((.., source_start..source_end))
            .map_err(h5ad_error)?;
        if values.shape() != [rows, source_end - source_start] {
            return Err(Error::InvalidAtlas(
                "legacy H5AD X slice has an unexpected shape".into(),
            ));
        }
        for (row, source_row) in values.outer_iter().enumerate() {
            let output_row = &mut output[row];
            for &(output_column, source_column) in &selected_columns {
                output_row[output_column] = source_row[source_column - source_start];
            }
        }
    }
    if output.iter().flatten().any(|value| !value.is_finite()) {
        return Err(Error::InvalidAtlas(
            "legacy H5AD X contains non-finite expression values".into(),
        ));
    }
    Ok(output)
}

fn legacy_column_groups(
    indices: &[usize],
    chunk_columns: usize,
) -> BTreeMap<usize, Vec<(usize, usize)>> {
    let mut groups: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
    for (output_column, &source_column) in indices.iter().enumerate() {
        groups
            .entry(source_column / chunk_columns)
            .or_default()
            .push((output_column, source_column));
    }
    groups
}

fn required_string_column(
    frame: &DataFrame,
    candidates: &[&str],
    description: &str,
) -> Result<Vec<String>> {
    optional_string_column(frame, candidates).ok_or_else(|| {
        let available = frame
            .get_column_names()
            .iter()
            .map(|name| name.as_str())
            .collect::<Vec<_>>()
            .join(", ");
        Error::InvalidAtlas(format!(
            "H5AD obs has no recognized {description} column; available columns: {available}"
        ))
    })
}

fn optional_string_column(frame: &DataFrame, candidates: &[&str]) -> Option<Vec<String>> {
    let column_name = candidates.iter().find_map(|candidate| {
        let normalized = normalize_column(candidate);
        frame
            .get_column_names()
            .into_iter()
            .find(|name| normalize_column(name.as_str()) == normalized)
            .map(|name| name.as_str().to_owned())
    })?;
    let column = frame.column(&column_name).ok()?;
    let strings = column.cast(&DataType::String).ok()?;
    let strings = strings.str().ok()?;
    Some(
        (0..frame.height())
            .map(|index| strings.get(index).unwrap_or("").to_owned())
            .collect(),
    )
}

fn normalize_column(name: &str) -> String {
    name.chars()
        .filter(|character| character.is_ascii_alphanumeric())
        .flat_map(char::to_lowercase)
        .collect()
}

fn matrix_to_rows(matrix: ArrayData, rows: usize, columns: usize) -> Result<Vec<Vec<f32>>> {
    let mut output = vec![vec![0.0; columns]; rows];
    match matrix {
        ArrayData::Array(matrix) => fill_dense(matrix, &mut output, rows, columns)?,
        ArrayData::CsrMatrix(matrix) => fill_csr(matrix, &mut output)?,
        ArrayData::CsrNonCanonical(matrix) => {
            let canonical = matrix.canonicalize().map_err(|_| {
                Error::InvalidAtlas("H5AD X contains an invalid sparse CSR matrix".into())
            })?;
            fill_csr(canonical, &mut output)?;
        }
        ArrayData::CscMatrix(matrix) => fill_csc(matrix, &mut output)?,
        ArrayData::DataFrame(_) => {
            return Err(Error::InvalidAtlas(
                "H5AD X is a data frame rather than a numeric matrix".into(),
            ));
        }
    }
    if output.iter().flatten().any(|value| !value.is_finite()) {
        return Err(Error::InvalidAtlas(
            "H5AD X contains non-finite expression values".into(),
        ));
    }
    Ok(output)
}

macro_rules! fill_dense_numeric {
    ($matrix:expr, $output:expr, $rows:expr, $columns:expr) => {{
        if $matrix.ndim() != 2 || $matrix.shape() != [$rows, $columns] {
            return Err(Error::InvalidAtlas(
                "H5AD X slice has an unexpected shape".into(),
            ));
        }
        for (offset, value) in $matrix.iter().enumerate() {
            $output[offset / $columns][offset % $columns] = *value as f32;
        }
    }};
}

fn fill_dense(
    matrix: DynArray,
    output: &mut [Vec<f32>],
    rows: usize,
    columns: usize,
) -> Result<()> {
    match matrix {
        DynArray::I8(matrix) => fill_dense_numeric!(matrix, output, rows, columns),
        DynArray::I16(matrix) => fill_dense_numeric!(matrix, output, rows, columns),
        DynArray::I32(matrix) => fill_dense_numeric!(matrix, output, rows, columns),
        DynArray::I64(matrix) => fill_dense_numeric!(matrix, output, rows, columns),
        DynArray::U8(matrix) => fill_dense_numeric!(matrix, output, rows, columns),
        DynArray::U16(matrix) => fill_dense_numeric!(matrix, output, rows, columns),
        DynArray::U32(matrix) => fill_dense_numeric!(matrix, output, rows, columns),
        DynArray::U64(matrix) => fill_dense_numeric!(matrix, output, rows, columns),
        DynArray::F32(matrix) => fill_dense_numeric!(matrix, output, rows, columns),
        DynArray::F64(matrix) => fill_dense_numeric!(matrix, output, rows, columns),
        DynArray::Bool(_) | DynArray::String(_) => {
            return Err(Error::InvalidAtlas("H5AD X is not numeric".into()));
        }
    }
    Ok(())
}

macro_rules! fill_csr_numeric {
    ($matrix:expr, $output:expr) => {{
        for (row_index, row) in $matrix.row_iter().enumerate() {
            for (&column_index, value) in row.col_indices().iter().zip(row.values()) {
                $output[row_index][column_index] = *value as f32;
            }
        }
    }};
}

fn fill_csr(matrix: DynCsrMatrix, output: &mut [Vec<f32>]) -> Result<()> {
    match matrix {
        DynCsrMatrix::I8(matrix) => fill_csr_numeric!(matrix, output),
        DynCsrMatrix::I16(matrix) => fill_csr_numeric!(matrix, output),
        DynCsrMatrix::I32(matrix) => fill_csr_numeric!(matrix, output),
        DynCsrMatrix::I64(matrix) => fill_csr_numeric!(matrix, output),
        DynCsrMatrix::U8(matrix) => fill_csr_numeric!(matrix, output),
        DynCsrMatrix::U16(matrix) => fill_csr_numeric!(matrix, output),
        DynCsrMatrix::U32(matrix) => fill_csr_numeric!(matrix, output),
        DynCsrMatrix::U64(matrix) => fill_csr_numeric!(matrix, output),
        DynCsrMatrix::F32(matrix) => fill_csr_numeric!(matrix, output),
        DynCsrMatrix::F64(matrix) => fill_csr_numeric!(matrix, output),
        DynCsrMatrix::Bool(_) | DynCsrMatrix::String(_) => {
            return Err(Error::InvalidAtlas("H5AD X is not numeric".into()));
        }
    }
    Ok(())
}

macro_rules! fill_csc_numeric {
    ($matrix:expr, $output:expr) => {{
        for (column_index, column) in $matrix.col_iter().enumerate() {
            for (&row_index, value) in column.row_indices().iter().zip(column.values()) {
                $output[row_index][column_index] = *value as f32;
            }
        }
    }};
}

fn fill_csc(matrix: DynCscMatrix, output: &mut [Vec<f32>]) -> Result<()> {
    match matrix {
        DynCscMatrix::I8(matrix) => fill_csc_numeric!(matrix, output),
        DynCscMatrix::I16(matrix) => fill_csc_numeric!(matrix, output),
        DynCscMatrix::I32(matrix) => fill_csc_numeric!(matrix, output),
        DynCscMatrix::I64(matrix) => fill_csc_numeric!(matrix, output),
        DynCscMatrix::U8(matrix) => fill_csc_numeric!(matrix, output),
        DynCscMatrix::U16(matrix) => fill_csc_numeric!(matrix, output),
        DynCscMatrix::U32(matrix) => fill_csc_numeric!(matrix, output),
        DynCscMatrix::U64(matrix) => fill_csc_numeric!(matrix, output),
        DynCscMatrix::F32(matrix) => fill_csc_numeric!(matrix, output),
        DynCscMatrix::F64(matrix) => fill_csc_numeric!(matrix, output),
        DynCscMatrix::Bool(_) | DynCscMatrix::String(_) => {
            return Err(Error::InvalidAtlas("H5AD X is not numeric".into()));
        }
    }
    Ok(())
}

fn h5ad_error(error: impl std::fmt::Display) -> Error {
    Error::InvalidAtlas(error.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legacy_columns_are_grouped_by_physical_chunk() {
        let groups = legacy_column_groups(&[108, 4, 110, 106], 107);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[&0], vec![(1, 4), (3, 106)]);
        assert_eq!(groups[&1], vec![(0, 108), (2, 110)]);
    }
}
