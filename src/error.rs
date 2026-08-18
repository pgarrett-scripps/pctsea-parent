use std::fmt::{Display, Formatter};
use std::io;
use std::path::PathBuf;

/// Errors produced by atlas I/O, input validation, or analysis.
#[derive(Debug)]
pub enum Error {
    Io(io::Error),
    InvalidAtlas(String),
    InvalidInput {
        path: PathBuf,
        line: usize,
        message: String,
    },
    InvalidConfig(String),
    Network(String),
    UnknownArgument(String),
}

pub type Result<T> = std::result::Result<T, Error>;

impl Display for Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(error) => write!(f, "I/O error: {error}"),
            Self::InvalidAtlas(message) => write!(f, "invalid atlas: {message}"),
            Self::InvalidInput {
                path,
                line,
                message,
            } => write!(f, "{}:{line}: {message}", path.display()),
            Self::InvalidConfig(message) => write!(f, "invalid configuration: {message}"),
            Self::Network(message) => write!(f, "network error: {message}"),
            Self::UnknownArgument(message) => write!(f, "{message}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            _ => None,
        }
    }
}

impl From<io::Error> for Error {
    fn from(value: io::Error) -> Self {
        Self::Io(value)
    }
}
