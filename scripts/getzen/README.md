# getzen

A small Rust utility for downloading all files from a [Zenodo](https://zenodo.org) record given its record ID. Handles concurrent downloads, checksum verification, resume-on-restart, and automatic retries on transient failures.

## Features

- Downloads every file attached to a Zenodo record in one command
- Concurrent downloads with a configurable parallelism limit
- Streams files to disk (safe for multi-GB records)
- Verifies MD5/SHA-256 checksums after download
- Skips files that already exist with a matching checksum (resume support)
- Automatic retry on network or API failures (up to 3 attempts with backoff)
- Supports restricted records via a Zenodo personal access token
- Works against both production Zenodo and Zenodo Sandbox

## Requirements

- Rust 1.85 or newer (edition 2024 is required by transitive dependencies)
- Install via [rustup](https://rustup.rs) if not already available:
  ```bash
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
  ```

## Build

```bash
cargo build --release
```

The compiled binary lands at `./target/release/getzen`. You can copy it anywhere on your `PATH`, or install it system-wide with:

```bash
cargo install --path .
```

## Usage

```bash
getzen <RECORD_ID> [OPTIONS]
```

### Examples

Download all files from a record into `./1234567/`:
```bash
getzen 1234567
```

Custom output directory with 8 concurrent downloads:
```bash
getzen 1234567 --output ./data --concurrency 8
```

Access a restricted record using a personal access token:
```bash
ZENODO_TOKEN=your_token_here getzen 1234567
```

Use the Zenodo Sandbox environment (for testing):
```bash
getzen 1234567 --sandbox
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-o, --output <DIR>` | Output directory | `./<record_id>` |
| `-c, --concurrency <N>` | Number of parallel downloads | `4` |
| `-r, --resume` | Skip files already present with matching checksum | `true` |
| `-t, --token <TOKEN>` | Zenodo access token (or set `ZENODO_TOKEN`) | — |
| `--sandbox` | Use `sandbox.zenodo.org` instead of production | `false` |

## How to find a record ID

Every Zenodo record has a URL like `https://zenodo.org/records/1234567`. The trailing number is the record ID.

## Retry behavior

Both the metadata request and each file download are wrapped in an automatic retry loop:

| Attempt | Wait before attempt |
|---------|--------------------:|
| 1       | none                |
| 2       | 1 second            |
| 3       | 5 seconds           |

After the third failure the operation gives up and reports the last error. Retries are per-file, so one file failing does not block the others. Failed attempts print a `[retry]` line to stderr so you can see what's happening.

Note: retries restart a file from byte 0 rather than resuming mid-stream. For flaky connections with very large files, running the tool again with `--resume` (on by default) will skip files already completed successfully.

## Notes

- Zenodo's public API is rate-limited (100 requests per minute per IP for anonymous requests, higher for authenticated). Keep `--concurrency` modest for large records.
- Checksum verification uses whatever algorithm Zenodo reports (`md5:` or `sha256:` prefix). Unknown algorithms produce a warning rather than a hard failure.
- Interrupted downloads: files with a matching checksum on disk are skipped on the next run, but partially written files will be re-downloaded from scratch (no HTTP `Range` support yet).

## License

MIT
