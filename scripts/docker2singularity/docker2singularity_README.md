# Docker to Singularity Converter

A Python script that converts Docker `run` commands to Singularity `exec` commands, automatically mapping Docker flags to their Singularity counterparts.

## Features

- **Volume mounts**: Converts `-v` or `--volume` to `--bind` or `-B`
- **Environment variables**: Converts `-e` or `--env` to `--env`
- **Working directory**: Converts `-w` or `--workdir` to `--pwd`
- **Image URI conversion**: Automatically adds `docker://` prefix for Docker Hub images
- **Multi-line formatting**: Optional formatting with backslash continuation for readability
- **Handles complex commands**: Properly tokenizes commands with quotes and escapes

## Installation

The script is located in `bin/docker2singularity.py` and requires Python 3.6+.

Make it executable:
```bash
chmod +x bin/docker2singularity.py
```

## Usage

### Basic Usage

Convert a Docker command and print to stdout:
```bash
python bin/docker2singularity.py input_file.txt
```

Save output to a file:
```bash
python bin/docker2singularity.py input_file.txt output_file.txt
```

### Multi-line Formatting

Format the output as a multi-line command with backslash continuation:
```bash
python bin/docker2singularity.py input_file.txt --multiline
```

Specify maximum line length (default is 80):
```bash
python bin/docker2singularity.py input_file.txt --multiline --max-line-length 100
```

### Help

View all options:
```bash
python bin/docker2singularity.py --help
```

## Examples

### Example 1: Basic Volume Mount

**Input file (docker_cmd.txt):**
```bash
docker run -v $PWD:/opt2 snakemake/snakemake:v5.24.2 \
    /opt2/chrom-seek run --assay cfChIP --genome hg19 --input \
    /opt2/.tests/sample.fastq.gz --output /opt2/output --dry-run
```

**Command:**
```bash
python bin/docker2singularity.py docker_cmd.txt --multiline
```

**Output:**
```bash
singularity exec --bind $PWD:/opt2 docker://snakemake/snakemake:v5.24.2 \
        /opt2/chrom-seek run --assay cfChIP --genome hg19 --input \
        /opt2/.tests/sample.fastq.gz --output /opt2/output --dry-run
```

### Example 2: With Environment Variables

**Input file:**
```bash
docker run -v /data:/mnt/data -e USER=admin -e PASSWORD=secret -w /workspace ubuntu:20.04 /bin/bash
```

**Output:**
```bash
singularity exec --bind /data:/mnt/data --env USER=admin --env PASSWORD=secret --pwd /workspace docker://ubuntu:20.04 /bin/bash
```

### Example 3: Multiple Volume Mounts

**Input file:**
```bash
docker run -v /data1:/mnt/data1 -v /data2:/mnt/data2 -v /config:/etc/config alpine:latest ls -la
```

**Output:**
```bash
singularity exec --bind /data1:/mnt/data1 --bind /data2:/mnt/data2 --bind /config:/etc/config docker://alpine:latest ls -la
```

## Supported Docker Flags

| Docker Flag | Singularity Equivalent | Status |
|-------------|------------------------|--------|
| `-v`, `--volume` | `--bind`, `-B` | ✅ Supported |
| `-e`, `--env` | `--env` | ✅ Supported |
| `-w`, `--workdir` | `--pwd` | ✅ Supported |
| `-u`, `--user` | N/A | ⚠️ Parsed but not mapped (Singularity runs as user by default) |
| `-it`, `-i`, `-t` | N/A | ⚠️ Ignored (interactive flags) |
| `--rm` | N/A | ⚠️ Ignored (not applicable) |
| `-d`, `--detach` | N/A | ⚠️ Ignored (not applicable) |

## Notes

- The script automatically converts Docker Hub image names (e.g., `ubuntu:20.04`) to Docker URIs (`docker://ubuntu:20.04`)
- Line continuations (`\`) in the input are handled correctly
- Quoted arguments are preserved
- The script focuses on `docker run` → `singularity exec` conversion for executing commands, not for building containers

## Limitations

- Does not support all Docker flags (see table above)
- Does not handle Docker-compose files
- Does not convert `docker build` commands
- Port mappings (`-p`) are not supported (Singularity handles networking differently)
- Some advanced Docker features may not have direct Singularity equivalents

## Requirements

- Python 3.6 or higher
- No external dependencies (uses only standard library)

## Contributing

If you encounter issues or need support for additional flags, please submit an issue or pull request.
