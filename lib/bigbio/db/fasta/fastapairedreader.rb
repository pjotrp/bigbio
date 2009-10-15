# FASTA paired reader keeps track of two FASTA files containing
# matching NT and AA sequences.
#

class FastaPairedReader

  def initialize ntfn, aafn, opts={:regex => '(\S+)'}
    @nt = FastaReader.new(ntfn, opts)
    @aa = FastaReader.new(aafn, opts)
  end

  # return a NT+AA pair
  def get id
    nt = @nt.get(id)
    aa = @aa.get(id)
    FastaPairedRecord.new(nt, aa)
  end

end
