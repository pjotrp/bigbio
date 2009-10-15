# Paired FASTA writer (tracks matching NT and AA sequences in two
# FASTA files)
#

class FastaPairedWriter

  def initialize ntfn, aafn
    @nt = FastaWriter.new(ntfn)
    @aa = FastaWriter.new(aafn)
  end

  def write rec
    @nt.write rec.nt
    @aa.write rec.aa
  end

  def close
    @nt.close
    @aa.close
  end
end
