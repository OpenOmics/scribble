# Code-snippets


Remove a files contents while keeping it's modification time:

```bash
function erase(){
    # Deletes a file and creates
    # a new file with the same
    # access and modification 
    # time of the deleted file
    # @INPUT: file(s)
    # @Returns: None
    function parse(){
        # Parses output of stat
        # command to get atime,
        # mtime, and ctime
        # @INPUT: $1=file
        # @INPUT: $2=pattern
        # @RETURN: parsed time    
        stat "${1}" \
            | grep -i --color=never "^${2}:" \
            | tail -n1 \
            | awk -F  '.' '{print $1}' \
            | awk '{print $(NF-1)$(NF)}' \
            | sed 's/-//g; s/://1; s/:/./1;'
    }
    
    for f in "$@"; do
        # Check for correct usage,
        # requires file(s) as input
        [ -d "${f}" ] \
            && echo "Warning: skipping directory... ${f}" \
            && continue
        [ ! -f "${f}" ] \
            && echo "Warning: skipping non-existent file... ${f}" \
            && continue
        atime=$(parse "${f}" "Access")
        mtime=$(parse "${f}" "Modify")
        rm "${f}"
        # Create empty file
        touch  "${f}"
        # Update access and 
        # modification times
        touch -a -t "${atime}" "${f}"
        touch -m -t "${mtime}" "${f}";
    done

}

export erase
```

**Usage**: `erase /path/to/files/*.txt`
