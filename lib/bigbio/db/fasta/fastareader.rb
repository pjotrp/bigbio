# Indexed FastaReader
#

require 'bigbio/db/fasta/indexer'

class FastaReader

  include Indexer

  # Initalize the reader of FASTA file _fn_. Options can be :regex and
  # :index (true/false)
  def initialize fn, opts = {}
    @f = File.open(fn)
    @fread_once = false
    @regex = opts[:regex]
    @regex = '^(\S+)' if @regex == nil
    indexer_use opts[:index]
  end

  # Parse the FASTA file and yield id, descr, sequence. When the indexer is on
  # it will index the records the first time. Note that, with indexing, when
  # you don't complete parsing there will be an error the second time. This is
  # a  # trade-off, otherwise one would always have to index the file and read
  # it twice.
  def parse_each
    @f.seek 0    # force file rewind
    @rec_fpos = 0
    @rec_line = @f.gets
    fpos = 0
    @count = 0
    begin
      # digest id from record description
      id, descr = digest_tag(@rec_line)
      id_fpos = @rec_fpos
      # parse the sequence
      seq = ""
      begin
        fpos = @f.tell
        line = @f.gets
        break if line =~ /^>/
        seq += line.strip 
      end while !@f.eof 
      # new record
      @count += 1
      @rec_fpos = fpos
      @rec_line = line
      # p [@rec_line, id, id_fpos]
      indexer_set(id, id_fpos) if @indexer and not @fread_once
      yield id, descr, seq
    end while !@f.eof
    @fread_once = true
  end

  # returns a FastaRecord for every item (invokes parse_each)
  def each
    parse_each { | id, descr, seq | yield FastaRecord.new(id, descr, seq) }
  end

  def first
    parse_each { | id, descr, seq | 
      return FastaRecord.new(id, descr, seq) 
    }
  end

  # Return a record by its +id+, nil when not found
  def get id
    indexed?
    if fpos = indexer_get(id)
      get_rec(fpos)
    else
      nil
    end
  end

  def get_rec fpos
    @f.seek fpos
    tag = @f.gets
    seq = ""
    begin
      line = @f.gets
      break if line =~ /^>/
      seq += line.strip 
    end while !@f.eof
    id, descr = digest_tag(tag)
    FastaRecord.new(id,descr,seq)
  end

  def get_by_index idx
    indexed?
    if fpos = indexer_get_by_index(idx)[1]
      ret = get_rec(fpos)
      return ret
    end
    nil
  end

  def digest_tag tag
    if tag =~ /^>/
      descr = $'.strip
      if descr =~ /#{@regex}/
        id = $1
        # p [descr,id]
        return id, descr
      end
      p descr  # do not remove these
      p @regex
    end
    raise "Can not digest '#{tag}' using '"+@regex+"'"
  end

  # Returns the size of the dataset - as read. After the final
  # record the size represents the number of items in the FASTA file
  def size
    @count
  end

  def close
    @f.close
  end

  private

  def indexed?
    if @indexer and not @fread_once
      # force indexer
      # $stderr.print "Force indexer"
      parse_each { | x, y, z | nil }
    end
    true
  end

end
