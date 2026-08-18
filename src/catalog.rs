use crate::{Error, Result};
use std::env;
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

const HDF5_MAGIC: &[u8; 8] = b"\x89HDF\r\n\x1a\n";
const ZIP_MAGIC: &[u8; 4] = b"PK\x03\x04";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AtlasFileKind {
    H5ad,
    Zip,
}

/// A known, reproducible public atlas release.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AtlasRelease {
    pub name: &'static str,
    pub version: &'static str,
    pub file_name: &'static str,
    pub download_url: &'static str,
    pub landing_page: &'static str,
    pub doi: &'static str,
    pub license: &'static str,
    pub expected_bytes: u64,
    pub kind: AtlasFileKind,
}

/// Official Human Cell Landscape release used for the first H5AD backend.
pub const HUMAN_CELL_LANDSCAPE: AtlasRelease = AtlasRelease {
    name: "hcl",
    version: "figshare-v4",
    file_name: "HCL_Fig1_adata.h5ad",
    download_url: "https://ndownloader.figshare.com/files/17727365",
    landing_page: "https://figshare.com/articles/dataset/HCL_DGE_Data/7235471",
    doi: "10.6084/m9.figshare.7235471",
    license: "CC BY 4.0",
    expected_bytes: 830_846_460,
    kind: AtlasFileKind::H5ad,
};

/// Official cell-info workbook paired with [`HUMAN_CELL_LANDSCAPE`].
pub const HUMAN_CELL_LANDSCAPE_ANNOTATIONS: AtlasRelease = AtlasRelease {
    name: "hcl",
    version: "figshare-v4",
    file_name: "HCL_Fig1_cell_Info.xlsx",
    download_url: "https://ndownloader.figshare.com/files/21758835",
    landing_page: "https://figshare.com/articles/dataset/HCL_DGE_Data/7235471",
    doi: "10.6084/m9.figshare.7235471",
    license: "CC BY 4.0",
    expected_bytes: 19_772_723,
    kind: AtlasFileKind::Zip,
};

/// Result of ensuring that a release exists locally.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtlasDownload {
    pub path: PathBuf,
    pub bytes: u64,
    pub downloaded: bool,
    pub release: AtlasRelease,
}

/// Returns the default local path without creating directories.
///
/// `PCTSEA_DATA_DIR` takes priority, followed by the platform data directory.
/// The final fallback is `.pctsea` below the current working directory.
pub fn default_atlas_path(release: AtlasRelease) -> Result<PathBuf> {
    let root = if let Some(path) = env::var_os("PCTSEA_DATA_DIR") {
        PathBuf::from(path)
    } else if let Some(path) = env::var_os("XDG_DATA_HOME") {
        PathBuf::from(path).join("pctsea")
    } else if let Some(path) = env::var_os("LOCALAPPDATA") {
        PathBuf::from(path).join("pctsea")
    } else {
        env::current_dir()?.join(".pctsea")
    };
    Ok(root
        .join("atlases")
        .join(release.name)
        .join(release.version)
        .join(release.file_name))
}

/// Downloads a known release with resumable `.part` files and atomic completion.
///
/// The callback receives `(downloaded_bytes, expected_total_bytes)` after each
/// chunk. Existing completed files are reused unless `force` is true.
pub fn download_atlas<F>(
    release: AtlasRelease,
    output: impl AsRef<Path>,
    force: bool,
    mut progress: F,
) -> Result<AtlasDownload>
where
    F: FnMut(u64, Option<u64>),
{
    let output = output.as_ref();
    if output.is_file() && !force {
        validate_download(output, release)?;
        return Ok(AtlasDownload {
            path: output.to_path_buf(),
            bytes: output.metadata()?.len(),
            downloaded: false,
            release,
        });
    }
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    let part = part_path(output);
    if force {
        if output.exists() {
            fs::remove_file(output)?;
        }
        if part.exists() {
            fs::remove_file(&part)?;
        }
    }
    let resume_at = part.metadata().map_or(0, |metadata| metadata.len());
    let agent = ureq::Agent::new_with_defaults();
    let mut request = agent.get(release.download_url);
    if resume_at > 0 {
        request = request.header("Range", &format!("bytes={resume_at}-"));
    }
    let response = request
        .call()
        .map_err(|error| Error::Network(error.to_string()))?;
    let status = response.status().as_u16();
    if status != 200 && status != 206 {
        return Err(Error::Network(format!(
            "atlas server returned unexpected HTTP status {status}"
        )));
    }
    if status == 206 && resume_at == 0 {
        return Err(Error::Network(
            "atlas server returned a partial response without a resume request".into(),
        ));
    }
    let resumed = resume_at > 0 && status == 206;
    let starting_bytes = if resumed { resume_at } else { 0 };
    let remaining = response.body().content_length();
    let expected_total = remaining.and_then(|bytes| bytes.checked_add(starting_bytes));
    let mut output_file = if resumed {
        OpenOptions::new().create(true).append(true).open(&part)?
    } else {
        File::create(&part)?
    };
    let (_, body) = response.into_parts();
    let mut reader = body.into_reader();
    let mut buffer = vec![0_u8; 1024 * 1024];
    let mut downloaded = starting_bytes;
    progress(downloaded, expected_total);
    loop {
        let count = reader
            .read(&mut buffer)
            .map_err(|error| Error::Network(error.to_string()))?;
        if count == 0 {
            break;
        }
        output_file.write_all(&buffer[..count])?;
        downloaded += count as u64;
        progress(downloaded, expected_total);
    }
    output_file.flush()?;
    output_file.sync_all()?;
    if expected_total.is_some_and(|expected| downloaded != expected) {
        return Err(Error::Network(format!(
            "incomplete download: received {downloaded} of {} bytes; rerun to resume",
            expected_total.unwrap()
        )));
    }
    validate_download(&part, release)?;
    fs::rename(&part, output)?;
    write_provenance(output, release, downloaded)?;
    Ok(AtlasDownload {
        path: output.to_path_buf(),
        bytes: downloaded,
        downloaded: true,
        release,
    })
}

fn validate_download(path: &Path, release: AtlasRelease) -> Result<()> {
    let length = path.metadata()?.len();
    if length != release.expected_bytes {
        return Err(Error::InvalidAtlas(format!(
            "{} has {length} bytes, expected {} for {} {}",
            path.display(),
            release.expected_bytes,
            release.name,
            release.version,
        )));
    }
    let mut file = File::open(path)?;
    let mut magic = [0_u8; HDF5_MAGIC.len()];
    file.read_exact(&mut magic)?;
    let valid = match release.kind {
        AtlasFileKind::H5ad => &magic == HDF5_MAGIC,
        AtlasFileKind::Zip => &magic[..ZIP_MAGIC.len()] == ZIP_MAGIC,
    };
    if !valid {
        return Err(Error::InvalidAtlas(format!(
            "{} does not have the expected {:?} signature; the server may have returned an error page",
            path.display(),
            release.kind,
        )));
    }
    Ok(())
}

fn part_path(output: &Path) -> PathBuf {
    let mut name = output.as_os_str().to_os_string();
    name.push(".part");
    PathBuf::from(name)
}

fn write_provenance(output: &Path, release: AtlasRelease, bytes: u64) -> Result<()> {
    let mut name = output.as_os_str().to_os_string();
    name.push(".provenance.txt");
    let path = PathBuf::from(name);
    let contents = format!(
        "atlas={}\nversion={}\nfile={}\nbytes={}\ndownload_url={}\nlanding_page={}\ndoi={}\nlicense={}\n",
        release.name,
        release.version,
        release.file_name,
        bytes,
        release.download_url,
        release.landing_page,
        release.doi,
        release.license,
    );
    fs::write(path, contents)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_path_includes_version_and_filename() {
        let path = default_atlas_path(HUMAN_CELL_LANDSCAPE).unwrap();
        assert!(path.ends_with("hcl/figshare-v4/HCL_Fig1_adata.h5ad"));
    }

    #[test]
    fn part_file_is_adjacent() {
        assert_eq!(
            part_path(Path::new("/tmp/atlas.h5ad")),
            PathBuf::from("/tmp/atlas.h5ad.part")
        );
    }
}
