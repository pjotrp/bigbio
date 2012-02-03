# Fasta writer

class FastaWriter

  # Open a FASTA stream for writing
  def initialize fn
    @f = File.open(fn,"w")
  end

  # write a FASTA item. Normally write using id,seq. If only one item is
  # passed in, the item should respond to descr and seq, or id and seq. 
  def write id, seq = nil
    if seq == nil
      item = id
      if item.respond_to?(:descr)
        @f.write ">"+item.descr+"\n"
      else
        @f.write ">"+item.id+"\n"
      end
      @f.write item.seq.to_s.strip+"\n"
    else
      @f.write ">"+id+"\n"
      @f.write seq.to_s.strip+"\n"
    end
  end

  def close
    @f.close
  end

end
