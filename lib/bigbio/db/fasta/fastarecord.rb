
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
    if nt.seq.size == aa.seq.size*3-1
      # account for EMBOSS cleverness
      nt.seq.chop!
      nt.seq.chop!
      aa.seq.chop!
    end
    if nt.seq.size == aa.seq.size*3-2
      # account for EMBOSS cleverness
      nt.seq.chop!
      aa.seq.chop!
    end
    if nt.seq.size == aa.seq.size*3-3
      aa.seq.chop!
    end
    nt_size = nt.seq.size
    expected_size = aa.seq.size*3
    # raise "Sequence size mismatch for #{nt.id} <nt:#{nt.seq.size} != #{aa.seq.size*3} (aa:#{aa.seq.size}*3)>" if expected_size - 3 > nt_size and  nt_size > expected_size + 3
  end

  def id
    @aa.id
  end
end
