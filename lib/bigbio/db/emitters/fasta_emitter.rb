
module Bio
  module Big
    class FastaEmitter

      def initialize fn, max_size = 100000
        @fn = fn
        @max_size = max_size
      end

      # Yield sequence information in sections of a maximum
      # size - usually iterators load the full sequence, but
      # without penalty it is possible to use a lot less 
      # memory.
      def emit_seq 
        f = File.open(@fn)
        tag = tag_digest(f.gets.strip)
        seq = ""
        index = 0
        begin
          line = f.gets.strip
          if line[0] =~ /^>/
            yield :tail,index,tag,seq
            tag = tag_digest(line)
            seq = ""
            index += 1
          else
            seq += line
          end
          while seq.size > @max_size
            yield :mid,index,tag,seq[0..@max_size-1]
            seq = seq[@max_size..-1]
          end
        end while !f.eof
        yield :tail,index,tag,seq
      end

      def tag_digest tag 
        if tag[0] == '>'
          tag[1..-1]
        else
          raise "Tag error #{tag}"
        end
      end
    end

  end # Big
end # Bio
