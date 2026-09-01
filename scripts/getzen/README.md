# zenodo_dl

A small Rust utility for downloading all files from a [Zenodo](https://zenodo.org) record given its record ID. Handles concurrent downloads, checksum verification, and resume-on-restart.

## Features

- Downloads every file attached to a Zenodo record in one command
- Concurrent downloads with a configurable parallelism limit
- Streams files to disk (safe for multi-GB records)
- Verifies MD5/SHA-256 checksums after download
- Skips files that already exist with a matching checksum (resume support)
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

The compiled binary lands at `./target/release/zenodo_dl`. You can copy it anywhere on your `PATH`, or install it system-wide with:

```bash
cargo install --path .
```

## Usage

```bash
zenodo_dl <RECORD_ID> [OPTIONS]
```

### Examples

Download all files from a record into `./1234567/`:
```bash
zenodo_dl 1234567
```

Custom output directory with 8 concurrent downloads:
```bash
zenodo_dl 1234567 --output ./data --concurrency 8
```

Access a restricted record using a personal access token:
```bash
ZENODO_TOKEN=your_token_here zenodo_dl 1234567
```

Use the Zenodo Sandbox environment (for testing):
```bash
zenodo_dl 1234567 --sandbox
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

## Notes

- Zenodo's public API is rate-limited (100 requests per minute per IP for anonymous requests, higher for authenticated). Keep `--concurrency` modest for large records.
- Checksum verification uses whatever algorithm Zenodo reports (`md5:` or `sha256:` prefix). Unknown algorithms produce a warning rather than a hard failure.
- Interrupted downloads: files with a matching checksum on disk are skipped on the next run, but partially written files will be re-downloaded from scratch (no HTTP `Range` support yet).

## License

MIT (or whatever you prefer).
