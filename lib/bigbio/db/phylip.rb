# Simple phylip reader. Supports PAML style files formatted as
#
# sequence 1
# AAGCTTCACCGGCGCAGTCATTCTCATAAT
# CGCCCACGGACTTACATCCTCATTACTATT
# sequence 2
# AAGCTTCACCGGCGCAATTATCCTCATAAT
# CGCCCACGGACTTACATCCTCATTATTATT
# sequence 3
# AAGCTTCACCGGCGCAGTTGTTCTTATAAT
# TGCCCACGGACTTACATCATCATTATTATT
# sequence 4
# AAGCTTCACCGGCGCAACCACCCTCATGAT
# TGCCCATGGACTCACATCCTCCCTACTGTT

module Bio
  module Big
    module PhylipReader
      # Define get_line as a lambda function, e.g.
      #   Bio::Big::PhylipReader.emit_seq(-> { lines.next }) { | name, seq | p [name,seq] }

      def PhylipReader::emit_seq get_line
        line = get_line.call.strip
        a = line.split
        seq_num = a[0].to_i
        seq_size = a[1].to_i
        name = nil
        seq = ""
        while true
          line = get_line.call
          break if line == nil or line == ""
          line = line.strip
          if name == nil
            name = line
            next
          end
          seq += line
          if seq.size >= seq_size
            raise "Name wrong size for #{name}" if name.size > 20
            raise "Sequence wrong size for #{name}" if seq.size > seq_size
            yield name, seq
            name = nil
            seq = ""
          end
        end
      end
    end
  end
end
