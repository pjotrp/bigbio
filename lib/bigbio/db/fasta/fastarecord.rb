
class FastaRecord
  attr_accessor :id, :descr, :seq

  def initialize id, descr, seq
    @id = id
    @descr = descr
    @seq = seq
  end
end

class FastaPairedRecord
  attr_reader :nt, :aa

  def initialize nt, aa
    @nt = nt
    @aa = aa
    raise "ID error NT #{nt.id} not matching AA #{aa.id}" if nt.id != aa.id
    raise "Sequence size mismatch for #{nt.id}" if nt.seq.size != aa.seq.size*3
  end

  def id
    @aa.id
  end
end
