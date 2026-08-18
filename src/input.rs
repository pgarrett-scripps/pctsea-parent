use crate::{Error, Result};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// One gene and its quantitative value in the experimental query.
#[derive(Clone, Debug, PartialEq)]
pub struct QueryGene {
    pub name: String,
    pub value: f64,
}

/// A quantitative gene/protein query.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneQuery {
    pub genes: Vec<QueryGene>,
}

impl GeneQuery {
    /// Creates a query, normalizing gene identifiers to uppercase.
    pub fn new<I, S>(genes: I) -> Result<Self>
    where
        I: IntoIterator<Item = (S, f64)>,
        S: Into<String>,
    {
        let mut seen = HashSet::new();
        let mut result = Vec::new();
        for (name, value) in genes {
            let name = name.into().trim().to_uppercase();
            if name.is_empty() {
                return Err(Error::InvalidConfig(
                    "query contains an empty gene name".into(),
                ));
            }
            if !value.is_finite() {
                return Err(Error::InvalidConfig(format!(
                    "query value for {name} is not finite"
                )));
            }
            if !seen.insert(name.clone()) {
                return Err(Error::InvalidConfig(format!(
                    "query contains duplicate gene {name}"
                )));
            }
            result.push(QueryGene { name, value });
        }
        if result.is_empty() {
            return Err(Error::InvalidConfig("query is empty".into()));
        }
        Ok(Self { genes: result })
    }

    /// Reads `gene<TAB>value`. A first header row is optional.
    pub fn from_tsv(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let reader = BufReader::new(File::open(path)?);
        let mut values = Vec::new();
        for (index, line) in reader.lines().enumerate() {
            let line_number = index + 1;
            let line = line?;
            if line.trim().is_empty() || line.trim_start().starts_with('#') {
                continue;
            }
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() < 2 {
                return Err(invalid(path, line_number, "expected gene<TAB>value"));
            }
            let value = match fields[1].trim().parse::<f64>() {
                Ok(value) => value,
                Err(_) if values.is_empty() && is_header(fields[0], fields[1]) => continue,
                Err(_) => return Err(invalid(path, line_number, "value is not a number")),
            };
            values.push((fields[0].to_string(), value));
        }
        Self::new(values)
    }
}

fn is_header(first: &str, second: &str) -> bool {
    matches!(
        first.trim().to_ascii_lowercase().as_str(),
        "gene" | "protein" | "id"
    ) || matches!(
        second.trim().to_ascii_lowercase().as_str(),
        "value" | "expression" | "abundance" | "score"
    )
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

    #[test]
    fn rejects_duplicate_case_insensitively() {
        let result = GeneQuery::new([("actb", 1.0), ("ACTB", 2.0)]);
        assert!(matches!(result, Err(Error::InvalidConfig(_))));
    }
}
