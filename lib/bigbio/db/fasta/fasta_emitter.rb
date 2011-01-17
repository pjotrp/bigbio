
module Bio
  module Big
    class FastaEmitter

      def initialize fn
        @fn = fn
      end

      # Yield sequence information in sections of a maximum
      # size
      def emit_seq max_size=10000
        f = File.open(@fn)
        tag = tag_digest(f.gets.strip)
        seq = ""
        begin
          line = f.gets.strip
          if line[0] =~ /^>/
            yield :tail,tag,seq
            tag = tag_digest(line)
            seq = ""
          else
            seq += line
          end
          while seq.size > max_size
            yield :mid,tag,seq[0..max_size-1]
            seq = seq[max_size..-1]
          end
        end while !f.eof
        yield :tail,tag,seq
      end

      def digest_tag tag 
        if tag.shift == '>'
          tag[1..-1]
        else
          raise "Tag error #{tag}"
        end
      end
    end

  end # Big
end # Bio
