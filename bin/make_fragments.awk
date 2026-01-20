# Fragment a FASTA into fixed windows.
# Params: -v LEN=50 -v STEP=25 -v PREFIX="DB"
# Output headers are safe-ish, but you should still sanitize before cons/subfam.
BEGIN{
  FS=""
  OFS=""
  if(LEN=="") LEN=50
  if(STEP=="") STEP=25
  if(PREFIX=="") PREFIX="DB"
  name=""; seq=""
}
function flush(    L,i,s,e,frag,hn){
  if(name=="") return
  gsub(/[ \t\r\n]/,"",seq)
  L=length(seq)
  if(L < LEN){ seq=""; return }
  for(i=1; i<=L-LEN+1; i+=STEP){
    s=i
    e=i+LEN-1
    frag=substr(seq,s,LEN)
    hn=name
    gsub(/[ \t\r\n]/,"_",hn)
    gsub(/[^A-Za-z0-9_.-]/,"_",hn)
    printf(">%s|%s|%d-%d\n%s\n", PREFIX, hn, s, e, frag)
  }
  seq=""
}
(/^>/){
  flush()
  name=substr($0,2)
  seq=""
  next
}
{
  line=$0
  gsub(/[ \t\r\n]/,"",line)
  seq=seq line
}
END{ flush() }
