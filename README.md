# FMI SDP course project: Huffman code compression
## Features
- Compress and archive files and directories
- \* and ? supported as wildcards to select multiple files and directories
- List files in an archive
- Check an archive or individual files in an archive for integrity
- Update single file in archive (faster than recreating the archive)
- Uses [Canonical Huffman Code](https://en.wikipedia.org/wiki/Canonical_Huffman_code) for more space-efficient storage

## Usage
huffman-coder [OPTION]... [FILE]...
- -f, --file  
specify the archive file
- -a, --add-file  
add file, ignoring special symbols, useful for files beginning with dash, or containing '*' and '?'
- -C, --directory  
change directory before operation
- -r, --raw  
only compress file content, not saving directory metadata
- -o, --output  
raw mode output file (default stdout)
- -c, --create  
create an archive from files
- -x, --extract  
extract files from an archive
- -l, --list  
list files in an archive
- -u, --update  
update a file in an archive
- -i, --check-integirty  
check an archive for integrity, and lists corrupted files
- -s, --add-integrity  
adds integrity data to the archive
- -e, --add-integrity-each  
adds integrity data to the archive, and for each file
- -L, --add-compression-levels  
adds compression levels to each file
- -U, --save-frequencies  
save frequencies to the archive to make updating faster