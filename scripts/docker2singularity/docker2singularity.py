#!/usr/bin/env python3
"""
Convert Docker run commands to Singularity exec commands.

This script reads a text file containing a Docker run command and converts it
to the equivalent Singularity exec command, mapping Docker flags to their
Singularity counterparts.

Usage:
    python docker2singularity.py <input_file> [output_file]
    
If output_file is not specified, the converted command is printed to stdout.
"""

import sys
import re
import argparse
from typing import List


class DockerToSingularityConverter:
    """Converts Docker run commands to Singularity exec commands."""
    
    def __init__(self):
        self.volume_mounts = []
        self.environment_vars = []
        self.working_dir = None
        self.user = None
        self.image = None
        self.command_args = []
        self.other_flags = []
        
    def parse_docker_command(self, command_text: str) -> bool:
        """
        Parse a Docker run command and extract its components.
        
        Args:
            command_text: The Docker run command as a string
            
        Returns:
            True if parsing was successful, False otherwise
        """
        # Remove line continuations and extra whitespace
        command_text = re.sub(r'\\\s*\n\s*', ' ', command_text)
        command_text = ' '.join(command_text.split())
        
        # Check if this is a docker run command
        if not re.search(r'\bdocker\s+run\b', command_text):
            print("Error: Not a valid 'docker run' command", file=sys.stderr)
            return False
        
        # Extract the command after 'docker run'
        match = re.search(r'docker\s+run\s+(.*)', command_text)
        if not match:
            return False
            
        args_string = match.group(1)
        
        # Tokenize the command, respecting quotes
        tokens = self._tokenize(args_string)
        
        i = 0
        while i < len(tokens):
            token = tokens[i]
            
            # Volume mounts: -v or --volume
            if token in ['-v', '--volume']:
                if i + 1 < len(tokens):
                    self.volume_mounts.append(tokens[i + 1])
                    i += 2
                else:
                    i += 1
                    
            # Environment variables: -e or --env
            elif token in ['-e', '--env']:
                if i + 1 < len(tokens):
                    self.environment_vars.append(tokens[i + 1])
                    i += 2
                else:
                    i += 1
                    
            # Working directory: -w or --workdir
            elif token in ['-w', '--workdir']:
                if i + 1 < len(tokens):
                    self.working_dir = tokens[i + 1]
                    i += 2
                else:
                    i += 1
                    
            # User: -u or --user
            elif token in ['-u', '--user']:
                if i + 1 < len(tokens):
                    self.user = tokens[i + 1]
                    i += 2
                else:
                    i += 1
                    
            # Flags that don't have direct Singularity equivalents but are common
            elif token in ['-it', '-i', '-t', '--rm', '-d', '--detach']:
                # Skip these flags as they don't apply to Singularity exec
                i += 1
                
            # Image name (first argument without a dash that's not after a flag)
            elif not token.startswith('-') and self.image is None:
                self.image = token
                # Everything after the image is the command to execute
                self.command_args = tokens[i + 1:]
                break
                
            # Other flags (capture but may not use)
            else:
                if token.startswith('-'):
                    self.other_flags.append(token)
                i += 1
                
        return self.image is not None
    
    def _tokenize(self, text: str) -> List[str]:
        """
        Tokenize command line arguments, respecting quotes.
        
        Args:
            text: Command line string to tokenize
            
        Returns:
            List of tokens
        """
        tokens = []
        current_token = []
        in_quotes = None
        escape_next = False
        
        for char in text:
            if escape_next:
                current_token.append(char)
                escape_next = False
                continue
                
            if char == '\\':
                escape_next = True
                continue
                
            if char in ['"', "'"]:
                if in_quotes == char:
                    in_quotes = None
                elif in_quotes is None:
                    in_quotes = char
                else:
                    current_token.append(char)
                continue
                
            if char.isspace() and in_quotes is None:
                if current_token:
                    tokens.append(''.join(current_token))
                    current_token = []
                continue
                
            current_token.append(char)
            
        if current_token:
            tokens.append(''.join(current_token))
            
        return tokens
    
    def convert_to_singularity(self) -> str:
        """
        Convert the parsed Docker command to Singularity exec command.
        
        Returns:
            The Singularity exec command as a string
        """
        if not self.image:
            return ""
        
        singularity_parts = ["singularity exec"]
        
        # Add volume mounts (--bind or -B)
        for mount in self.volume_mounts:
            singularity_parts.append(f"--bind {mount}")
        
        # Add environment variables (--env)
        for env_var in self.environment_vars:
            singularity_parts.append(f"--env {env_var}")
        
        # Add working directory (--pwd)
        if self.working_dir:
            singularity_parts.append(f"--pwd {self.working_dir}")
        
        # Convert image name to docker:// URI if it doesn't have a scheme
        image = self.image
        if not re.match(r'^\w+://', image):
            # It's a Docker Hub image
            image = f"docker://{image}"
        
        singularity_parts.append(image)
        
        # Add the command and its arguments
        if self.command_args:
            singularity_parts.extend(self.command_args)
        
        return ' '.join(singularity_parts)
    
    def format_multiline(self, command: str, max_line_length: int = 80) -> str:
        """
        Format a long command into multiple lines with backslash continuation.
        
        Args:
            command: The command string to format
            max_line_length: Maximum length of each line
            
        Returns:
            Multi-line formatted command
        """
        parts = command.split()
        lines = []
        current_line = []
        current_length = 0
        
        for part in parts:
            part_length = len(part) + 1  # +1 for space
            
            if current_length + part_length > max_line_length and current_line:
                lines.append(' '.join(current_line) + ' \\')
                current_line = [part]
                current_length = len(part)
            else:
                current_line.append(part)
                current_length += part_length
        
        if current_line:
            lines.append(' '.join(current_line))
        
        # Indent continuation lines
        formatted_lines = [lines[0]]
        for line in lines[1:]:
            formatted_lines.append('        ' + line)
        
        return '\n'.join(formatted_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Convert Docker run commands to Singularity exec commands",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python docker2singularity.py docker_command.txt
    python docker2singularity.py docker_command.txt singularity_command.txt
    python docker2singularity.py input.txt --multiline
        """
    )
    
    parser.add_argument(
        'input_file',
        type=str,
        help='Input file containing Docker run command'
    )
    
    parser.add_argument(
        'output_file',
        type=str,
        nargs='?',
        help='Output file for Singularity exec command (optional, defaults to stdout)'
    )
    
    parser.add_argument(
        '--multiline',
        action='store_true',
        help='Format output as multi-line command with backslash continuation'
    )
    
    parser.add_argument(
        '--max-line-length',
        type=int,
        default=80,
        help='Maximum line length for multiline formatting (default: 80)'
    )
    
    args = parser.parse_args()
    
    # Read the Docker command from input file
    try:
        with open(args.input_file, 'r') as f:
            docker_command = f.read()
    except FileNotFoundError:
        print(f"Error: File '{args.input_file}' not found", file=sys.stderr)
        return 1
    except IOError as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        return 1
    
    # Convert the command
    converter = DockerToSingularityConverter()
    
    if not converter.parse_docker_command(docker_command):
        print("Error: Failed to parse Docker command", file=sys.stderr)
        return 1
    
    singularity_command = converter.convert_to_singularity()
    
    if not singularity_command:
        print("Error: Failed to convert command", file=sys.stderr)
        return 1
    
    # Format if requested
    if args.multiline:
        singularity_command = converter.format_multiline(
            singularity_command, 
            args.max_line_length
        )
    
    # Output the result
    if args.output_file:
        try:
            with open(args.output_file, 'w') as f:
                f.write(singularity_command)
                f.write('\n')
            print(f"Converted command written to '{args.output_file}'")
        except IOError as e:
            print(f"Error writing file: {e}", file=sys.stderr)
            return 1
    else:
        print(singularity_command)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
