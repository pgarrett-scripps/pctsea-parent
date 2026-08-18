use crate::{Error, GeneQuery, Result};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

const MAGIC: &[u8; 8] = b"PCTSEA\0\x01";
const FORMAT_VERSION: u32 = 1;
const MAX_STRING_BYTES: usize = 16 * 1024 * 1024;
const MAX_RECORDS: usize = 1_000_000_000;

/// Metadata for one cell in an atlas.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Cell {
    pub name: String,
    pub cell_type: String,
    pub dataset: String,
}

/// A non-zero expression value returned from a gene lookup.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneExpression {
    pub cell_index: usize,
    pub value: f32,
}

/// Summary of a native atlas.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtlasInfo {
    pub format: String,
    pub format_version: Option<u32>,
    pub cells: usize,
    pub genes: usize,
    pub non_zero_expressions: Option<u64>,
    pub cell_types: usize,
    pub datasets: usize,
}

#[derive(Clone, Debug)]
struct GeneIndex {
    name: String,
    offset: u64,
    entries: u32,
}

/// An opened native atlas. Only metadata and the gene index are resident.
#[derive(Debug)]
pub struct Atlas {
    path: PathBuf,
    cells: Vec<Cell>,
    genes: Vec<GeneIndex>,
    gene_lookup: HashMap<String, usize>,
}

/// Query-aligned expression values loaded from an atlas.
#[derive(Clone, Debug)]
pub struct LoadedQuery {
    pub gene_names: Vec<String>,
    pub query_values: Vec<f64>,
    /// `cell_values[cell_index][gene_index]`.
    pub cell_values: Vec<Vec<f32>>,
    pub missing_genes: Vec<String>,
}

/// Backend-neutral access required by the analysis engine.
///
/// Native and H5AD atlases implement the same contract, keeping storage
/// decisions out of the statistics.
pub trait AtlasData {
    fn cells(&self) -> &[Cell];
    fn info(&self) -> AtlasInfo;
    fn gene_expression(&self, gene: &str) -> Result<Option<Vec<GeneExpression>>>;
    fn load_query(&self, query: &GeneQuery) -> Result<LoadedQuery>;
}

impl Atlas {
    /// Opens and validates a native `.pctsea` atlas.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let mut reader = BufReader::new(File::open(&path)?);
        let mut magic = [0_u8; 8];
        reader.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(Error::InvalidAtlas("bad magic bytes".into()));
        }
        let version = read_u32(&mut reader)?;
        if version != FORMAT_VERSION {
            return Err(Error::InvalidAtlas(format!(
                "unsupported format version {version}"
            )));
        }
        let cell_count = checked_count(read_u32(&mut reader)?, "cell")?;
        let gene_count = checked_count(read_u32(&mut reader)?, "gene")?;
        let mut cells = Vec::with_capacity(cell_count);
        for _ in 0..cell_count {
            cells.push(Cell {
                name: read_string(&mut reader)?,
                cell_type: read_string(&mut reader)?,
                dataset: read_string(&mut reader)?,
            });
        }
        let mut genes = Vec::with_capacity(gene_count);
        let mut gene_lookup = HashMap::with_capacity(gene_count);
        for index in 0..gene_count {
            let name = read_string(&mut reader)?;
            let offset = read_u64(&mut reader)?;
            let entries = read_u32(&mut reader)?;
            if entries as usize > cell_count {
                return Err(Error::InvalidAtlas(format!(
                    "gene {name} has more postings than cells"
                )));
            }
            if gene_lookup.insert(name.clone(), index).is_some() {
                return Err(Error::InvalidAtlas(format!("duplicate gene {name}")));
            }
            genes.push(GeneIndex {
                name,
                offset,
                entries,
            });
        }
        let file_length = reader.get_ref().metadata()?.len();
        for gene in &genes {
            let end = gene
                .offset
                .checked_add(u64::from(gene.entries) * 8)
                .ok_or_else(|| Error::InvalidAtlas("posting offset overflow".into()))?;
            if end > file_length {
                return Err(Error::InvalidAtlas(format!(
                    "postings for {} extend past end of file",
                    gene.name
                )));
            }
        }
        Ok(Self {
            path,
            cells,
            genes,
            gene_lookup,
        })
    }

    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    pub fn contains_gene(&self, gene: &str) -> bool {
        self.gene_lookup.contains_key(&gene.trim().to_uppercase())
    }

    pub fn info(&self) -> AtlasInfo {
        AtlasInfo {
            format: "pctsea".into(),
            format_version: Some(FORMAT_VERSION),
            cells: self.cells.len(),
            genes: self.genes.len(),
            non_zero_expressions: Some(self.genes.iter().map(|gene| u64::from(gene.entries)).sum()),
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

    /// Reads the sparse postings for one gene without loading other expressions.
    pub fn gene_expression(&self, gene: &str) -> Result<Option<Vec<GeneExpression>>> {
        let key = gene.trim().to_uppercase();
        let Some(&index) = self.gene_lookup.get(&key) else {
            return Ok(None);
        };
        let gene = &self.genes[index];
        let mut reader = BufReader::new(File::open(&self.path)?);
        reader.seek(SeekFrom::Start(gene.offset))?;
        let mut values = Vec::with_capacity(gene.entries as usize);
        let mut previous = None;
        for _ in 0..gene.entries {
            let cell_index = read_u32(&mut reader)? as usize;
            let value = f32::from_le_bytes(read_array(&mut reader)?);
            if cell_index >= self.cells.len() {
                return Err(Error::InvalidAtlas(format!(
                    "gene {} refers to missing cell {cell_index}",
                    gene.name
                )));
            }
            if previous.is_some_and(|last| last >= cell_index) {
                return Err(Error::InvalidAtlas(format!(
                    "postings for {} are not strictly sorted",
                    gene.name
                )));
            }
            if !value.is_finite() || value == 0.0 {
                return Err(Error::InvalidAtlas(format!(
                    "gene {} contains an invalid sparse value",
                    gene.name
                )));
            }
            previous = Some(cell_index);
            values.push(GeneExpression { cell_index, value });
        }
        Ok(Some(values))
    }

    /// Loads only the expression columns present in `query`.
    pub fn load_query(&self, query: &GeneQuery) -> Result<LoadedQuery> {
        let mut gene_names = Vec::new();
        let mut query_values = Vec::new();
        let mut postings = Vec::new();
        let mut missing_genes = Vec::new();
        for query_gene in &query.genes {
            match self.gene_expression(&query_gene.name)? {
                Some(values) => {
                    gene_names.push(query_gene.name.clone());
                    query_values.push(query_gene.value);
                    postings.push(values);
                }
                None => missing_genes.push(query_gene.name.clone()),
            }
        }
        if gene_names.is_empty() {
            return Err(Error::InvalidConfig(
                "none of the query genes occur in the atlas".into(),
            ));
        }
        let mut cell_values = vec![vec![0_f32; gene_names.len()]; self.cells.len()];
        for (gene_index, values) in postings.into_iter().enumerate() {
            for expression in values {
                cell_values[expression.cell_index][gene_index] = expression.value;
            }
        }
        Ok(LoadedQuery {
            gene_names,
            query_values,
            cell_values,
            missing_genes,
        })
    }
}

impl AtlasData for Atlas {
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

/// Converts portable long-form TSV files into the indexed native format.
pub struct AtlasBuilder;

impl AtlasBuilder {
    /// Packs a cell table and an expression table.
    ///
    /// `cells.tsv`: `cell<TAB>cell_type<TAB>dataset` (`dataset` optional).
    /// `expressions.tsv`: `gene<TAB>cell<TAB>expression`.
    pub fn pack_tsv(
        cells_path: impl AsRef<Path>,
        expressions_path: impl AsRef<Path>,
        output_path: impl AsRef<Path>,
    ) -> Result<AtlasInfo> {
        let cells_path = cells_path.as_ref();
        let expressions_path = expressions_path.as_ref();
        let cells = read_cells(cells_path)?;
        if cells.len() > u32::MAX as usize {
            return Err(Error::InvalidConfig(
                "too many cells for atlas format".into(),
            ));
        }
        let lookup: HashMap<&str, u32> = cells
            .iter()
            .enumerate()
            .map(|(index, cell)| (cell.name.as_str(), index as u32))
            .collect();
        let expressions = read_expressions(expressions_path, &lookup)?;
        write_atlas(output_path.as_ref(), &cells, expressions)?;
        Atlas::open(output_path).map(|atlas| atlas.info())
    }
}

fn read_cells(path: &Path) -> Result<Vec<Cell>> {
    let reader = BufReader::new(File::open(path)?);
    let mut cells = Vec::new();
    let mut names = HashSet::new();
    for (index, line) in reader.lines().enumerate() {
        let number = index + 1;
        let line = line?;
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if cells.is_empty()
            && fields.len() >= 2
            && fields[0].trim().eq_ignore_ascii_case("cell")
            && fields[1].trim().to_ascii_lowercase().contains("type")
        {
            continue;
        }
        if fields.len() < 2 {
            return Err(invalid(
                path,
                number,
                "expected cell<TAB>cell_type[<TAB>dataset]",
            ));
        }
        let name = fields[0].trim().to_string();
        let cell_type = fields[1].trim().to_string();
        let dataset = fields.get(2).map_or("", |value| value.trim()).to_string();
        if name.is_empty() || cell_type.is_empty() {
            return Err(invalid(
                path,
                number,
                "cell and cell_type must not be empty",
            ));
        }
        if !names.insert(name.clone()) {
            return Err(invalid(path, number, "duplicate cell name"));
        }
        cells.push(Cell {
            name,
            cell_type,
            dataset,
        });
    }
    if cells.is_empty() {
        return Err(invalid(path, 0, "cell table is empty"));
    }
    Ok(cells)
}

fn read_expressions(
    path: &Path,
    cells: &HashMap<&str, u32>,
) -> Result<BTreeMap<String, Vec<(u32, f32)>>> {
    let reader = BufReader::new(File::open(path)?);
    let mut genes: BTreeMap<String, Vec<(u32, f32)>> = BTreeMap::new();
    for (index, line) in reader.lines().enumerate() {
        let number = index + 1;
        let line = line?;
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 3 {
            return Err(invalid(
                path,
                number,
                "expected gene<TAB>cell<TAB>expression",
            ));
        }
        if genes.is_empty()
            && fields[0].trim().eq_ignore_ascii_case("gene")
            && fields[1].trim().eq_ignore_ascii_case("cell")
        {
            continue;
        }
        let gene = fields[0].trim().to_uppercase();
        let cell = fields[1].trim();
        let Some(&cell_index) = cells.get(cell) else {
            return Err(invalid(path, number, &format!("unknown cell {cell}")));
        };
        let value = fields[2]
            .trim()
            .parse::<f32>()
            .map_err(|_| invalid(path, number, "expression is not a 32-bit number"))?;
        if gene.is_empty() || !value.is_finite() {
            return Err(invalid(
                path,
                number,
                "gene is empty or expression is not finite",
            ));
        }
        if value != 0.0 {
            genes.entry(gene).or_default().push((cell_index, value));
        }
    }
    if genes.is_empty() {
        return Err(invalid(path, 0, "expression table has no non-zero values"));
    }
    for postings in genes.values_mut() {
        postings.sort_unstable_by_key(|&(cell, _)| cell);
        let mut merged: Vec<(u32, f32)> = Vec::with_capacity(postings.len());
        for &(cell, value) in postings.iter() {
            if let Some((last_cell, last_value)) = merged.last_mut()
                && *last_cell == cell
            {
                *last_value += value;
            } else {
                merged.push((cell, value));
            }
        }
        if merged.iter().any(|(_, value)| !value.is_finite()) {
            return Err(invalid(
                path,
                0,
                "duplicate expressions overflowed while being summed",
            ));
        }
        merged.retain(|(_, value)| *value != 0.0);
        *postings = merged;
    }
    genes.retain(|_, postings| !postings.is_empty());
    Ok(genes)
}

fn write_atlas(
    path: &Path,
    cells: &[Cell],
    expressions: BTreeMap<String, Vec<(u32, f32)>>,
) -> Result<()> {
    if expressions.len() > u32::MAX as usize {
        return Err(Error::InvalidConfig(
            "too many genes for atlas format".into(),
        ));
    }
    let mut writer = BufWriter::new(File::create(path)?);
    writer.write_all(MAGIC)?;
    write_u32(&mut writer, FORMAT_VERSION)?;
    write_u32(&mut writer, cells.len() as u32)?;
    write_u32(&mut writer, expressions.len() as u32)?;
    for cell in cells {
        write_string(&mut writer, &cell.name)?;
        write_string(&mut writer, &cell.cell_type)?;
        write_string(&mut writer, &cell.dataset)?;
    }
    let mut offset_positions = Vec::with_capacity(expressions.len());
    for (gene, postings) in &expressions {
        write_string(&mut writer, gene)?;
        offset_positions.push(writer.stream_position()?);
        write_u64(&mut writer, 0)?;
        write_u32(&mut writer, postings.len() as u32)?;
    }
    let mut offsets = Vec::with_capacity(expressions.len());
    for postings in expressions.values() {
        offsets.push(writer.stream_position()?);
        for &(cell, value) in postings {
            write_u32(&mut writer, cell)?;
            writer.write_all(&value.to_le_bytes())?;
        }
    }
    for (position, offset) in offset_positions.into_iter().zip(offsets) {
        writer.seek(SeekFrom::Start(position))?;
        write_u64(&mut writer, offset)?;
    }
    writer.flush()?;
    Ok(())
}

fn checked_count(value: u32, name: &str) -> Result<usize> {
    let value = value as usize;
    if value > MAX_RECORDS {
        Err(Error::InvalidAtlas(format!("unreasonable {name} count")))
    } else {
        Ok(value)
    }
}

fn read_string(reader: &mut impl Read) -> Result<String> {
    let length = read_u32(reader)? as usize;
    if length > MAX_STRING_BYTES {
        return Err(Error::InvalidAtlas("unreasonable string length".into()));
    }
    let mut bytes = vec![0_u8; length];
    reader.read_exact(&mut bytes)?;
    String::from_utf8(bytes).map_err(|_| Error::InvalidAtlas("string is not UTF-8".into()))
}

fn write_string(writer: &mut impl Write, value: &str) -> Result<()> {
    if value.len() > u32::MAX as usize {
        return Err(Error::InvalidConfig("string is too long".into()));
    }
    write_u32(writer, value.len() as u32)?;
    writer.write_all(value.as_bytes())?;
    Ok(())
}

fn read_array<const N: usize>(reader: &mut impl Read) -> Result<[u8; N]> {
    let mut bytes = [0_u8; N];
    reader.read_exact(&mut bytes)?;
    Ok(bytes)
}

fn read_u32(reader: &mut impl Read) -> Result<u32> {
    Ok(u32::from_le_bytes(read_array(reader)?))
}

fn read_u64(reader: &mut impl Read) -> Result<u64> {
    Ok(u64::from_le_bytes(read_array(reader)?))
}

fn write_u32(writer: &mut impl Write, value: u32) -> Result<()> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn write_u64(writer: &mut impl Write, value: u64) -> Result<()> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn invalid(path: &Path, line: usize, message: &str) -> Error {
    Error::InvalidInput {
        path: path.to_path_buf(),
        line,
        message: message.into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_file(name: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("pctsea-{nonce}-{name}"))
    }

    #[test]
    fn packs_opens_and_loads_selected_genes() {
        let cells = temp_file("cells.tsv");
        let expressions = temp_file("expressions.tsv");
        let atlas_path = temp_file("atlas.pctsea");
        fs::write(
            &cells,
            "cell\tcell_type\tdataset\nc1\tT cell\td1\nc2\tB cell\td1\n",
        )
        .unwrap();
        fs::write(
            &expressions,
            "gene\tcell\texpression\nACTB\tc1\t2\nactb\tc2\t3\nCD3D\tc1\t4\n",
        )
        .unwrap();
        let info = AtlasBuilder::pack_tsv(&cells, &expressions, &atlas_path).unwrap();
        assert_eq!(
            (info.cells, info.genes, info.non_zero_expressions),
            (2, 2, Some(3))
        );

        let atlas = Atlas::open(&atlas_path).unwrap();
        let query = GeneQuery::new([("actb", 10.0), ("missing", 2.0)]).unwrap();
        let loaded = atlas.load_query(&query).unwrap();
        assert_eq!(loaded.gene_names, ["ACTB"]);
        assert_eq!(loaded.cell_values, [vec![2.0], vec![3.0]]);
        assert_eq!(loaded.missing_genes, ["MISSING"]);

        for path in [cells, expressions, atlas_path] {
            fs::remove_file(path).unwrap();
        }
    }
}
