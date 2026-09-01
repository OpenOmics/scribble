// Cargo.toml dependencies:
// [dependencies]
// reqwest = { version = "0.12", features = ["stream", "json"] }
// tokio = { version = "1", features = ["full"] }
// serde = { version = "1", features = ["derive"] }
// serde_json = "1"
// futures-util = "0.3"
// indicatif = "0.17"
// clap = { version = "4", features = ["derive"] }
// anyhow = "1"
// sha2 = "0.10"
// md-5 = "0.10"

use anyhow::{anyhow, Context, Result};
use clap::Parser;
use futures_util::StreamExt;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use serde::Deserialize;
use std::future::Future;
use std::path::{Path, PathBuf};
use std::time::Duration;
use tokio::fs::File;
use tokio::io::AsyncWriteExt;

/// Delay (in seconds) to wait *before* each attempt.
/// Attempt 1: no wait. Attempt 2: 1s. Attempt 3: 5s.
const RETRY_DELAYS_SECS: &[u64] = &[0, 1, 5];

/// Retry a fallible async operation up to `RETRY_DELAYS_SECS.len()` times.
/// `label` is used for logging context on failed attempts.
async fn with_retries<T, E, F, Fut>(label: &str, mut op: F) -> Result<T, E>
where
    F: FnMut() -> Fut,
    Fut: Future<Output = Result<T, E>>,
    E: std::fmt::Display,
{
    let mut last_err: Option<E> = None;
    for (attempt, delay) in RETRY_DELAYS_SECS.iter().enumerate() {
        if *delay > 0 {
            tokio::time::sleep(Duration::from_secs(*delay)).await;
        }
        match op().await {
            Ok(v) => return Ok(v),
            Err(e) => {
                eprintln!(
                    "  [retry] {} — attempt {}/{} failed: {}",
                    label,
                    attempt + 1,
                    RETRY_DELAYS_SECS.len(),
                    e
                );
                last_err = Some(e);
            }
        }
    }
    Err(last_err.expect("loop always populates last_err on failure"))
}

#[derive(Parser, Debug)]
#[command(author, version, about = "Download all files from a Zenodo record")]
struct Args {
    /// Zenodo record ID (e.g. 1234567)
    record_id: String,

    /// Output directory (default: ./<record_id>)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Number of concurrent downloads
    #[arg(short, long, default_value_t = 4)]
    concurrency: usize,

    /// Skip files that already exist with matching checksum
    #[arg(short, long, default_value_t = true)]
    resume: bool,

    /// Optional Zenodo access token (for restricted records)
    #[arg(short, long, env = "ZENODO_TOKEN")]
    token: Option<String>,

    /// Use Zenodo Sandbox API instead of production
    #[arg(long, default_value_t = false)]
    sandbox: bool,
}

#[derive(Debug, Deserialize)]
struct ZenodoRecord {
    id: u64,
    #[serde(default)]
    title: String,
    files: Vec<ZenodoFile>,
}

#[derive(Debug, Deserialize, Clone)]
struct ZenodoFile {
    key: String,
    size: u64,
    checksum: String, // format: "md5:abc123..." (occasionally other algos)
    links: FileLinks,
}

#[derive(Debug, Deserialize, Clone)]
struct FileLinks {
    #[serde(rename = "self")]
    self_link: String,
}

async fn fetch_record(
    client: &reqwest::Client,
    record_id: &str,
    sandbox: bool,
    token: Option<&str>,
) -> Result<ZenodoRecord> {
    let base = if sandbox {
        "https://sandbox.zenodo.org"
    } else {
        "https://zenodo.org"
    };
    let url = format!("{}/api/records/{}", base, record_id);

    with_retries(&format!("fetch metadata for {}", record_id), || async {
        let mut req = client.get(&url);
        if let Some(t) = token {
            req = req.bearer_auth(t);
        }

        let resp = req.send().await.context("failed to reach Zenodo API")?;
        let status = resp.status();
        if !status.is_success() {
            let body = resp.text().await.unwrap_or_default();
            return Err(anyhow!(
                "Zenodo API returned {} for record {}: {}",
                status,
                record_id,
                body
            ));
        }

        let record: ZenodoRecord = resp.json().await.context("failed to parse record JSON")?;
        Ok(record)
    })
    .await
}

/// Verify a file on disk against a Zenodo checksum string like "md5:abcd..."
async fn verify_checksum(path: &Path, checksum: &str) -> Result<bool> {
    let (algo, expected) = checksum
        .split_once(':')
        .ok_or_else(|| anyhow!("malformed checksum: {}", checksum))?;

    let path = path.to_path_buf();
    let algo = algo.to_ascii_lowercase();
    let expected = expected.to_ascii_lowercase();

    // Hash in a blocking task since it's CPU-bound
    let actual = tokio::task::spawn_blocking(move || -> Result<String> {
        use std::io::Read;
        let mut file = std::fs::File::open(&path)?;
        let mut buf = [0u8; 64 * 1024];
        match algo.as_str() {
            "md5" => {
                use md5::{Digest, Md5};
                let mut hasher = Md5::new();
                loop {
                    let n = file.read(&mut buf)?;
                    if n == 0 { break; }
                    hasher.update(&buf[..n]);
                }
                Ok(format!("{:x}", hasher.finalize()))
            }
            "sha256" => {
                use sha2::{Digest, Sha256};
                let mut hasher = Sha256::new();
                loop {
                    let n = file.read(&mut buf)?;
                    if n == 0 { break; }
                    hasher.update(&buf[..n]);
                }
                Ok(format!("{:x}", hasher.finalize()))
            }
            other => Err(anyhow!("unsupported checksum algorithm: {}", other)),
        }
    })
    .await??;

    Ok(actual == expected)
}

async fn download_file(
    client: reqwest::Client,
    file: ZenodoFile,
    out_dir: PathBuf,
    token: Option<String>,
    resume: bool,
    multi: MultiProgress,
) -> Result<()> {
    let out_path = out_dir.join(&file.key);
    if let Some(parent) = out_path.parent() {
        tokio::fs::create_dir_all(parent).await.ok();
    }

    // Skip if already present with matching checksum
    if resume && out_path.exists() {
        if let Ok(true) = verify_checksum(&out_path, &file.checksum).await {
            let pb = multi.add(ProgressBar::new(file.size));
            pb.set_style(ProgressStyle::with_template("{msg}").unwrap());
            pb.finish_with_message(format!("✓ {} (cached)", file.key));
            return Ok(());
        }
    }

    let pb = multi.add(ProgressBar::new(file.size));
    pb.set_style(
        ProgressStyle::with_template(
            "{msg:<40} [{bar:30.cyan/blue}] {bytes}/{total_bytes} ({bytes_per_sec}, {eta})",
        )
        .unwrap()
        .progress_chars("=>-"),
    );
    pb.set_message(file.key.clone());

    // Retry the whole request+stream as a unit. On mid-stream failure we
    // truncate and start over rather than trying to resume with Range.
    let download_result = with_retries(&format!("download {}", file.key), || async {
        let mut req = client.get(&file.links.self_link);
        if let Some(t) = &token {
            req = req.bearer_auth(t);
        }

        let resp = req
            .send()
            .await
            .with_context(|| format!("request failed for {}", file.key))?;
        if !resp.status().is_success() {
            return Err(anyhow!(
                "download of {} failed: HTTP {}",
                file.key,
                resp.status()
            ));
        }

        // Reset progress + output file for this attempt
        pb.set_position(0);
        let mut out = File::create(&out_path)
            .await
            .with_context(|| format!("cannot create {}", out_path.display()))?;
        let mut stream = resp.bytes_stream();

        while let Some(chunk) = stream.next().await {
            let chunk = chunk.context("stream error")?;
            out.write_all(&chunk).await?;
            pb.inc(chunk.len() as u64);
        }
        out.flush().await?;
        Ok(())
    })
    .await;

    if let Err(e) = download_result {
        pb.finish_with_message(format!("✗ {} (failed)", file.key));
        return Err(e);
    }

    // Verify checksum after download
    match verify_checksum(&out_path, &file.checksum).await {
        Ok(true) => pb.finish_with_message(format!("✓ {}", file.key)),
        Ok(false) => {
            pb.finish_with_message(format!("✗ {} (checksum mismatch)", file.key));
            return Err(anyhow!("checksum mismatch for {}", file.key));
        }
        Err(e) => {
            pb.finish_with_message(format!("? {} (checksum error: {})", file.key, e));
        }
    }

    Ok(())
}

#[tokio::main]
async fn main() -> Result<()> {
    let args = Args::parse();

    let client = reqwest::Client::builder()
        .user_agent(concat!("zenodo-downloader/", env!("CARGO_PKG_VERSION")))
        .build()?;

    println!("Fetching metadata for record {}...", args.record_id);
    let record = fetch_record(&client, &args.record_id, args.sandbox, args.token.as_deref()).await?;

    println!("Record {}: \"{}\"", record.id, record.title);
    println!("{} file(s) to download\n", record.files.len());

    let out_dir = args.output.unwrap_or_else(|| PathBuf::from(&args.record_id));
    tokio::fs::create_dir_all(&out_dir).await?;

    let multi = MultiProgress::new();
    let token = args.token.clone();

    // Bounded concurrency using buffer_unordered
    let results: Vec<Result<()>> = futures_util::stream::iter(record.files.into_iter().map(|f| {
        let client = client.clone();
        let out_dir = out_dir.clone();
        let token = token.clone();
        let multi = multi.clone();
        async move { download_file(client, f, out_dir, token, args.resume, multi).await }
    }))
    .buffer_unordered(args.concurrency)
    .collect()
    .await;

    let mut failures = 0usize;
    for r in &results {
        if let Err(e) = r {
            eprintln!("Error: {:#}", e);
            failures += 1;
        }
    }

    if failures > 0 {
        Err(anyhow!("{} file(s) failed to download", failures))
    } else {
        println!("\nAll files downloaded to {}", out_dir.display());
        Ok(())
    }
}