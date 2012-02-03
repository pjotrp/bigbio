# Fasta writer

class FastaWriter

  # Open a FASTA stream for writing
  def initialize fn
    @f = File.open(fn,"w")
  end

  # write a FASTA item. An itex should respond to descr and seq,
  # or id and seq
  def write item
    if item.respond_to?(:descr)
      @f.write ">"+item.descr+"\n"
    else
      @f.write ">"+item.id+"\n"
    end
    @f.write item.seq.to_s.strip+"\n"
  end

  def close
    @f.close
  end

end
