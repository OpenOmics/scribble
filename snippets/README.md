# Snippets


The snippets directory contains [code-snippets](https://en.wikipedia.org/wiki/Snippet_(programming)). These are usually small stand-alone pieces of re-usable code. Snippets can act like templates for repeating certain acts/patterns/tasks. Did you figure out how to do something really cool in R or python that you want to share? Scribble it down here!

### Remove a files contents while keeping it's modification time:

*What?* This bash function is useful for deleting the contents of large files associated with a Snakemake pipeline while keeping the original modification times. This will prevent any DAG related errors when running or dry-running the pipeline again. It will also prevent certain rules from running. Snakemake looks at mtime to see if an input to a rule as been updated. 

*Why?* This command can be used to reduce the storage footprint of a pipeline's working directory while preserving current dry-run functionality. 

**Please keep in mind:** this bash function will delete any files that are provided as input. And as so, please take steps to ensure that any files that are provided to this function are already archived (if they are important) prior to running this command.

#### **Usage**: `erase /path/to/files/*.txt`

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

export -f erase
```
