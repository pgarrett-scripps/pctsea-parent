use pctsea::{
    AnalysisConfig, AnalysisResult, AnyAtlas, AtlasBuilder, Error, GeneQuery, HUMAN_CELL_LANDSCAPE,
    HUMAN_CELL_LANDSCAPE_ANNOTATIONS, Result, ScoringMethod, analyze, default_atlas_path,
    download_atlas,
};
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::process::ExitCode;

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("pctsea: {error}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<()> {
    let mut arguments = env::args().skip(1);
    let Some(command) = arguments.next() else {
        print_help();
        return Ok(());
    };
    let rest: Vec<_> = arguments.collect();
    if command == "help" || command == "--help" || command == "-h" {
        print_help();
        return Ok(());
    }
    if command == "version" || command == "--version" || command == "-V" {
        println!("pctsea {}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    if command == "atlas" {
        return atlas_command(&rest);
    }
    let options = Options::parse(&rest)?;
    if options.flags.contains("help") {
        print_help();
        return Ok(());
    }
    match command.as_str() {
        "pack" => pack(&options),
        "info" => info(&options),
        "gene" => gene(&options),
        "analyze" => run_analysis(&options),
        _ => Err(Error::UnknownArgument(format!(
            "unknown command {command:?}; run `pctsea help`"
        ))),
    }
}

fn atlas_command(arguments: &[String]) -> Result<()> {
    let Some(action) = arguments.first().map(String::as_str) else {
        return Err(Error::UnknownArgument(
            "missing atlas action; use list, path, or download".into(),
        ));
    };
    if action == "list" {
        if arguments.len() != 1 {
            return Err(Error::UnknownArgument(
                "`atlas list` does not accept additional arguments".into(),
            ));
        }
        let release = HUMAN_CELL_LANDSCAPE;
        println!("name\tversion\tformat\tlicense\tdoi");
        println!(
            "{}\t{}\th5ad+xlsx\t{}\t{}",
            release.name, release.version, release.license, release.doi
        );
        return Ok(());
    }
    let Some(name) = arguments.get(1).map(String::as_str) else {
        return Err(Error::UnknownArgument(format!(
            "`atlas {action}` needs an atlas name; currently supported: hcl"
        )));
    };
    let release = match name {
        "hcl" => HUMAN_CELL_LANDSCAPE,
        _ => {
            return Err(Error::UnknownArgument(format!(
                "unknown atlas {name:?}; currently supported: hcl"
            )));
        }
    };
    let options = Options::parse(&arguments[2..])?;
    match action {
        "path" => {
            options.ensure_known(&[])?;
            println!("h5ad\t{}", default_atlas_path(release)?.display());
            println!(
                "annotations\t{}",
                default_atlas_path(HUMAN_CELL_LANDSCAPE_ANNOTATIONS)?.display()
            );
            Ok(())
        }
        "download" => {
            options.ensure_known(&["output", "force"])?;
            let path = options
                .get("output")
                .map(Path::new)
                .map(Path::to_path_buf)
                .map_or_else(|| default_atlas_path(release), Ok)?;
            eprintln!(
                "acquiring {} {} from Figshare...",
                release.name, release.version
            );
            let mut last_report = 0_u64;
            let result = download_atlas(
                release,
                &path,
                options.flags.contains("force"),
                |done, total| {
                    if done == 0
                        || done.saturating_sub(last_report) >= 64 * 1024 * 1024
                        || total == Some(done)
                    {
                        match total {
                            Some(total) => eprintln!("downloaded {} / {}", size(done), size(total)),
                            None => eprintln!("downloaded {}", size(done)),
                        }
                        last_report = done;
                    }
                },
            )?;
            let annotations_path = path.with_file_name(HUMAN_CELL_LANDSCAPE_ANNOTATIONS.file_name);
            eprintln!("acquiring paired HCL cell annotations from Figshare...");
            last_report = 0;
            let annotations = download_atlas(
                HUMAN_CELL_LANDSCAPE_ANNOTATIONS,
                &annotations_path,
                options.flags.contains("force"),
                |done, total| {
                    if done == 0
                        || done.saturating_sub(last_report) >= 4 * 1024 * 1024
                        || total == Some(done)
                    {
                        match total {
                            Some(total) => eprintln!("downloaded {} / {}", size(done), size(total)),
                            None => eprintln!("downloaded {}", size(done)),
                        }
                        last_report = done;
                    }
                },
            )?;
            eprintln!("validating atlas bundle and preparing metadata cache...");
            let atlas = AnyAtlas::open(&path)?;
            let atlas_info = atlas.info();
            if result.downloaded {
                println!("downloaded\t{}", result.path.display());
            } else {
                println!("already_present\t{}", result.path.display());
            }
            println!("bytes\t{}", result.bytes);
            println!("annotations\t{}", annotations.path.display());
            println!("annotation_bytes\t{}", annotations.bytes);
            println!("cells\t{}", atlas_info.cells);
            println!("cell_types\t{}", atlas_info.cell_types);
            println!("genes\t{}", atlas_info.genes);
            println!("source\t{}", result.release.landing_page);
            Ok(())
        }
        _ => Err(Error::UnknownArgument(format!(
            "unknown atlas action {action:?}; use list, path, or download"
        ))),
    }
}

fn pack(options: &Options) -> Result<()> {
    options.ensure_known(&["cells", "expressions", "output"])?;
    let cells = options.required("cells")?;
    let expressions = options.required("expressions")?;
    let output = options.required("output")?;
    let summary = AtlasBuilder::pack_tsv(cells, expressions, output)?;
    println!(
        "packed {} cells, {} genes, and {} non-zero expressions into {}",
        summary.cells,
        summary.genes,
        summary.non_zero_expressions.unwrap_or(0),
        output
    );
    Ok(())
}

fn info(options: &Options) -> Result<()> {
    options.ensure_known(&["atlas"])?;
    let atlas = AnyAtlas::open(options.required("atlas")?)?;
    let summary = atlas.info();
    println!("format\t{}", summary.format);
    println!(
        "format_version\t{}",
        summary
            .format_version
            .map_or_else(|| "unknown".into(), |value| value.to_string())
    );
    println!("cells\t{}", summary.cells);
    println!("cell_types\t{}", summary.cell_types);
    println!("datasets\t{}", summary.datasets);
    println!("genes\t{}", summary.genes);
    println!(
        "non_zero_expressions\t{}",
        summary
            .non_zero_expressions
            .map_or_else(|| "unknown".into(), |value| value.to_string())
    );
    Ok(())
}

fn gene(options: &Options) -> Result<()> {
    options.ensure_known(&["atlas", "gene", "limit"])?;
    let atlas = AnyAtlas::open(options.required("atlas")?)?;
    let gene = options.required("gene")?;
    let limit = options.parse_or("limit", usize::MAX)?;
    let Some(expressions) = atlas.gene_expression(gene)? else {
        return Err(Error::InvalidConfig(format!(
            "gene {gene:?} is not in the atlas"
        )));
    };
    println!("gene\tcell\tcell_type\tdataset\texpression");
    for expression in expressions.into_iter().take(limit) {
        let cell = &atlas.cells()[expression.cell_index];
        println!(
            "{}\t{}\t{}\t{}\t{}",
            gene.to_uppercase(),
            clean(&cell.name),
            clean(&cell.cell_type),
            clean(&cell.dataset),
            expression.value
        );
    }
    Ok(())
}

fn run_analysis(options: &Options) -> Result<()> {
    options.ensure_known(&[
        "atlas",
        "input",
        "output",
        "scores-output",
        "method",
        "min-genes",
        "threshold",
        "permutations",
        "seed",
        "dataset",
        "no-threshold",
        "include-zeros",
    ])?;
    let atlas = AnyAtlas::open(options.required("atlas")?)?;
    let query = GeneQuery::from_tsv(options.required("input")?)?;
    let mut config = AnalysisConfig::default();
    if let Some(method) = options.get("method") {
        config.method = method.parse::<ScoringMethod>()?;
    }
    config.min_genes_per_cell = options.parse_or("min-genes", config.min_genes_per_cell)?;
    config.permutations = options.parse_or("permutations", config.permutations)?;
    config.random_seed = options.parse_or("seed", config.random_seed)?;
    if options.flags.contains("no-threshold") {
        config.score_threshold = None;
    } else if let Some(threshold) = options.get("threshold") {
        config.score_threshold = Some(parse_value("threshold", threshold)?);
    }
    config.positive_pairs_only = !options.flags.contains("include-zeros");
    config.datasets = options
        .values
        .get("dataset")
        .into_iter()
        .flatten()
        .flat_map(|value| value.split(','))
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(str::to_string)
        .collect();

    let result = analyze(&atlas, &query, &config)?;
    let output = options.get("output").unwrap_or("-");
    write_results(&result, output)?;
    if let Some(path) = options.get("scores-output") {
        write_scores(&result, path)?;
    }
    eprintln!(
        "matched {} genes ({} missing); scored {}/{} cells; {} passed",
        result.matched_genes.len(),
        result.missing_genes.len(),
        result.cells_scored,
        result.cells_considered,
        result.cells_passing
    );
    Ok(())
}

fn write_results(result: &AnalysisResult, path: &str) -> Result<()> {
    let mut writer: Box<dyn Write> = if path == "-" {
        Box::new(BufWriter::new(std::io::stdout().lock()))
    } else {
        Box::new(BufWriter::new(File::create(Path::new(path))?))
    };
    writeln!(
        writer,
        "cell_type\tcells_total\tcells_passing\thypergeometric_p\thypergeometric_fdr\tenrichment_score\tnormalized_enrichment_score\tpermutation_p\tpermutation_fdr\ttop_genes"
    )?;
    for result in &result.cell_types {
        let top_genes = result
            .top_genes
            .iter()
            .map(|(gene, count)| format!("{gene}[{count}]"))
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            clean(&result.cell_type),
            result.cells_total,
            result.cells_passing,
            number(result.hypergeometric_p_value),
            number(result.hypergeometric_fdr),
            number(result.enrichment_score),
            number(result.normalized_enrichment_score),
            number(result.permutation_p_value),
            number(result.permutation_fdr),
            top_genes
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_scores(result: &AnalysisResult, path: &str) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "rank\tcell\tcell_type\tdataset\tscore\tgenes_used\tpasses_threshold"
    )?;
    for (index, score) in result.cell_scores.iter().enumerate() {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            index + 1,
            clean(&score.cell_name),
            clean(&score.cell_type),
            clean(&score.dataset),
            score.score,
            score.genes_used,
            score.passes_threshold
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.8}")
    } else {
        "NA".into()
    }
}

fn clean(value: &str) -> String {
    value.replace(['\t', '\r', '\n'], " ")
}

#[derive(Debug)]
struct Options {
    values: HashMap<String, Vec<String>>,
    flags: HashSet<String>,
}

impl Options {
    fn parse(arguments: &[String]) -> Result<Self> {
        let mut values: HashMap<String, Vec<String>> = HashMap::new();
        let mut flags = HashSet::new();
        let flag_names = ["help", "no-threshold", "include-zeros", "force"];
        let mut index = 0;
        while index < arguments.len() {
            let argument = &arguments[index];
            let Some(name) = argument.strip_prefix("--") else {
                return Err(Error::UnknownArgument(format!(
                    "unexpected positional argument {argument:?}"
                )));
            };
            if flag_names.contains(&name) {
                flags.insert(name.to_string());
                index += 1;
                continue;
            }
            let Some(value) = arguments.get(index + 1) else {
                return Err(Error::UnknownArgument(format!("--{name} needs a value")));
            };
            if value.starts_with("--") {
                return Err(Error::UnknownArgument(format!("--{name} needs a value")));
            }
            values
                .entry(name.to_string())
                .or_default()
                .push(value.clone());
            index += 2;
        }
        Ok(Self { values, flags })
    }

    fn get(&self, name: &str) -> Option<&str> {
        self.values
            .get(name)
            .and_then(|values| values.last())
            .map(String::as_str)
    }

    fn required(&self, name: &str) -> Result<&str> {
        self.get(name)
            .ok_or_else(|| Error::UnknownArgument(format!("missing required --{name}")))
    }

    fn parse_or<T>(&self, name: &str, default: T) -> Result<T>
    where
        T: std::str::FromStr,
    {
        self.get(name)
            .map_or(Ok(default), |value| parse_value(name, value))
    }

    fn ensure_known(&self, known: &[&str]) -> Result<()> {
        for name in self.values.keys().chain(&self.flags) {
            if !known.contains(&name.as_str()) && name != "help" {
                return Err(Error::UnknownArgument(format!("unknown option --{name}")));
            }
        }
        Ok(())
    }
}

fn parse_value<T>(name: &str, value: &str) -> Result<T>
where
    T: std::str::FromStr,
{
    value
        .parse()
        .map_err(|_| Error::InvalidConfig(format!("--{name} has invalid value {value:?}")))
}

fn size(bytes: u64) -> String {
    const UNITS: &[&str] = &["B", "KiB", "MiB", "GiB", "TiB"];
    let mut value = bytes as f64;
    let mut unit = 0;
    while value >= 1024.0 && unit + 1 < UNITS.len() {
        value /= 1024.0;
        unit += 1;
    }
    format!("{value:.1} {}", UNITS[unit])
}

fn print_help() {
    println!(
        "pctsea - file-backed proteomic cell-type enrichment\n\n\
USAGE:\n  \
pctsea --version\n  \
pctsea atlas list\n  \
pctsea atlas path hcl\n  \
pctsea atlas download hcl [--output FILE] [--force]\n  \
pctsea pack --cells CELLS.tsv --expressions EXPRESSIONS.tsv --output ATLAS.pctsea\n  \
pctsea info --atlas ATLAS.h5ad\n  \
pctsea gene --atlas ATLAS.h5ad --gene GENE [--limit N]\n  \
pctsea analyze --atlas ATLAS.h5ad --input QUERY.tsv [OPTIONS]\n\n\
ANALYZE OPTIONS:\n  \
--output FILE             Cell-type results TSV; default is stdout\n  \
--scores-output FILE      Optional ranked cell scores TSV\n  \
--method METHOD           pearson (default), cosine, or dot-product\n  \
--min-genes N             Minimum matched genes per cell; default 3\n  \
--threshold X             Score threshold; default 0\n  \
--no-threshold            Retain every scorable cell\n  \
--permutations N          Label permutations; default 1000\n  \
--seed N                  Deterministic random seed\n  \
--dataset NAME            Restrict dataset; repeat or comma-separate\n  \
--include-zeros           Score zero/non-positive pairs too"
    );
}
