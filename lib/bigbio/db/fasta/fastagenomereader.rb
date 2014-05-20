# Buffered FastaGenomeReader
#

class BufferMissed < Exception
end

class FastaGenomeReader

  # Initalize the reader of FASTA file 
  def initialize fn, parse_descriptor_func, bufsize=64_000
    @f = File.open(fn)
    @parse_descriptor = parse_descriptor_func
    @bufsize = bufsize
    @buf = read_next
  end

  # Returns the reference nucleotide. When the buffer is missed a BufferMissed
  # exception is thrown.
  def ref chr,pos
  end

private

  # Fill the next buffer until the next descriptor is reached or the buffer
  # is full
  def read_next
    while (line = @f.gets)
      p line
    end
  end
end
