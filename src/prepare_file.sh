# Provide gff3 file as first argument
# File without a header
tail -n +2 $1 | grep gene > genes_only.gff3
# Replace naming column with ID only
cat genes_only.gff3 | sed -E 's/ID=([a-zA-Z0-9\.]*)(;Name=.*)?/\1/g' > genes_only_clean_ID.gff3