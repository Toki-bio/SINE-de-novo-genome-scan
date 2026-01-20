# Sanitize FASTA headers to safe ASCII: [A-Za-z0-9_.-]
# Usage: awk -f sanitize_fasta_headers.awk in.fa > out.fa
/^>/{
  h = substr($0,2)
  gsub(/[ \t\r\n]+/, "_", h)
  gsub(/[^A-Za-z0-9_.-]/, "_", h)
  gsub(/_+/, "_", h)
  if(h=="") h="EMPTY"
  print ">" h
  next
}
{ print }
