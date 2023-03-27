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

### Create a HPC datashare link

*What?* This bash function will create a shareable HPC datashare link. The datashare directory is a special directory which are readable via the web, but are not browseable. Here is more information from [HPC staff](https://hpc.nih.gov/nih/datashare.html). The datashare directory is not a graveyard, so please clean up after yourself if the data does not need to be shared. If you are moving files into the datashare directory, please DO NOT create duplicates of any data. Create hard links with the `ln` command to the original data. 

*Why?* This is useful for temporarily sharing data with collaborators or for creating shareable IGV sessions with other people. You can add this snippet to your `~/.bashrc` or `~/.bash_profile` for later use.

**Please keep in mind:** The datashare directory is not a permanent archival location. Please take steps to delete any data in this directory after you are done. The great thing about creating a link, is that now the data is downloadable with `wget` or `cURL`. If you are sharing BAM files, please make should their BAI files are in the same directory.

#### **Usage**: `hpclink /data/NCBR/datashare/files/*.fastq.gz`

```bash
function hpclink () {
    # Creates a HPC datashare link,
    # given a path to a file in the 
    # datashare directory, by default,
    # it is setup to create URLs pointing
    # /data/NCBR/datashare; however, this
    # can be edited to point to other 
    # datashare locations on Biowulf.
    # $@ = Data to share 
    # @RETURNS = Helix datashare links  
    local prefix="/gpfs/gsfs11/users/NCBR/datashare/"; 
    for f in "$@"; do 
        abspath="$(readlink -e "$f")"; 
        if [[ "$abspath" =~ ^"$prefix" ]]; then 
            link=$(
                echo "$abspath" \
                | sed "s@^$prefix@https:\/\/hpc.nih.gov/\~NCBR/@g"
            ); 
            echo "$link"; 
        fi; 
    done
}

# Export function 
export -f hpclink
```

### Embed an image into html

*What?* This bash function will take an image (such as a PNG) and embed it into html. One or more image can be provided. Each image is converted into base64 and the html is redirected to standard output.

*Why?* This is useful for quickly bundling up multiple images into a single file that can be opened with a modern browser. You can add this snippet to your `~/.bashrc` or `~/.bash_profile` for later use.

#### **Usage**: `image2html file1.png file2.png file3.png`

```bash
function image2html(){
    # Create html with an embeded image,
    # given one or more images.
    # $@ = One or more images
    # @RETURNS = Image embedded HTML
    for f in $@; do 
        b=$(cat "$f" | base64 -w 0); 
        echo "<img src=\"data:image/${f##*.};base64,${b}\" alt=\"${f%.*}\">"; 
    done
}

# Export fuction 
export -f image2html
```

### Extract files of various types

*What?* This bash function will extract/uncompress various file types. I can never remember the options for tar, so this is super useful.

*Why?* This [relevant xkcd](https://xkcd.com/1168/) explains everything. You can add this snippet to your `~/.bashrc` or `~/.bash_profile` for later use.

#### **Usage**: `extract archive.tar.gz`

```bash
extract () {
     # Extract files of various Types,
     # $1 = File to extract
     # @RETURNS: None, just extracts file
     if [ -f $1 ] ; then
         case $1 in
             *.tar.bz2)   tar vxjf $1    ;;
             *.tar.gz)    tar vxzf $1    ;;
             *.bz2)       bunzip2 $1     ;;
             *.rar)       rar x $1       ;;
             *.gz)        gunzip $1      ;;
             *.tar)       tar vxf $1     ;;
             *.tbz2)      tar vxjf $1    ;;
             *.tgz)       tar vxzf $1    ;;
             *.zip)       unzip $1       ;;
             *.Z)         uncompress $1  ;;
             *.7z)        7z x $1    ;;
             *)           echo "'$1' cannot be extracted via extract()" ;;
         esac
     else
         echo "'$1' is not a valid file"
     fi
}

# Export fuction 
export -f extract 
```


### Filter out scRNA-seq cells with outlier features

*What?* This snippet of R code will filter out cells containing outlier features-- such as nCount_RNA, nFeature_RNA and percent mitochondria-- in a list of seurat objects using median absolute deviations (MAD).

*Author* @OluwayioseOA

```R
# Filter out cells with outlier features,
# such as nCount_RNA, nFeature_RNA and % 
# mitochondria in a list of seurat objects 
# using median absolute deviations (MAD).

# Variable `sgcList` is a list of multiple
# seurat object SC data
for (i in 1:length(sgcList)) { 
    print(paste0("Sample ", i))
    sgc.sample <- sgcList[[i]]
    low_lib_size <- isOutlier(sgc.sample$nCount_RNA, log=TRUE, type="lower")
    high_lib_size <- isOutlier(sgc.sample$nCount_RNA, log=TRUE, type="higher")
    low_nfeature <- isOutlier(sgc.sample$nFeature_RNA, log=TRUE, type="lower")
    high_nfeature <- isOutlier(sgc.sample$nFeature_RNA, log=TRUE, type="higher")
    high_mit <- isOutlier(sgc.sample$percent.mt, log=TRUE, type="higher")
    sgc.sample <- sgc.sample %>% mutate(
        low_lib_size  = low_lib_size,
        high_lib_size = high_lib_size,
        low_nfeature  = low_nfeature,
        high_nfeature = high_nfeature,
        high_mit = high_mit
    )
    # Filter out the outliers
    sgc.sample<- subset(
        sgc.sample,
        low_lib_size == "FALSE" &
        high_lib_size =="FALSE" &
        low_nfeature == "FALSE" &
        high_nfeature == "FALSE" &
        high_mit == "FALSE"
    )
    sgcList[[i]] <- sgc.sample
}
```
