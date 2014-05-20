# Buffered FastaGenomeReader. 
#
# The logic is straightforward. Keep reading a file for the first matching
# chr,pos, max_bufsize characters at a time. Looking back is not allowed.
#

class FastaGenomeReader

  attr_accessor :rec

  class Record
    attr_reader :chr, :start, :stop, :descr
    attr_accessor :buf
    def initialize line =nil, func = nil
      @descr = line
      @chr,@start,@stop = func.call(line) if func
      @start = 0 if not @start
      @stop = -1 if not @stop
      @offset = @start
      @buf = ''
    end

    def value pos
      @buf[pos-@offset]
    end

    def in_range? chr,pos
      return false if chr != @chr
      @offset <= pos and pos < @offset+@buf.size
    end

    def move_offset
      @offset += @buf.size
    end

    def empty_buf
      move_offset
      @buf = ""
    end
  end

  # Initalize the reader of FASTA file 
  def initialize fn, parse_descriptor_func, max_bufsize=64_000
    @f = File.open(fn)
    @parse_descriptor_func = parse_descriptor_func
    @max_bufsize = max_bufsize
    @rec = Record.new
    read_next
  end

  # Returns the reference nucleotide. Chr can be any name (i.e., chr or a bin).
  def ref chr,pos
    while not @rec.in_range?(chr,pos)
      return nil if not read_next 
    end
    @rec.value(pos)
  end

  # Fetch in current chromosome/bin
  def [] pos
    ref(@rec.chr,pos)
  end

private

  # Fill the next buffer until the next descriptor is reached or the buffer
  # is full
  def read_next
    @rec.empty_buf
    while (line = @f.gets)
      next if line =~ /^#/
      if line =~ /^>/
        @rec = Record.new(line,@parse_descriptor_func)
      else
        @rec.buf << line.strip
        return @rec if @rec.buf.size > @max_bufsize
      end
    end
    nil
  end
end
