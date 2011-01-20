
require 'set'

module Bio
  module Big

    STOP_CODONS = Set.new(%w{TAG TAA TGA UAG UAA UGA})
    START_CODONS = Set.new(%w{ATG AUG})

    # The short frame uses a simpler concept. The sequence is immutable,
    # always forward and in frame 0. That makes it easy to reason. It 
    # also return all ORF's in one go, with the left/right locations.
    class ShortFrameState
      def initialize seq, min_size = 30
        @seq = seq.upcase  
        @min_size = (min_size/3).to_i+1
        @codons = @seq.scan(/(\w\w\w)/).flatten
      end

      def get_stopstop_orfs 
        list = get_codon_orfs1(Proc.new { | codon | STOP_CODONS.include?(codon) })
        list.map { |codons| 
          codons[1..-1].join
        }
      end

      def get_startstop_orfs 
        list = get_codon_orfs2(
                 Proc.new { | codon | START_CODONS.include?(codon) },
                 Proc.new { | codon | STOP_CODONS.include?(codon) })
        p list
        list.map { |codons| 
          codons.join
        }
      end

      def get_codon_orfs1 is_splitter_func, checkfirst=true
        orfs = split(@codons,is_splitter_func)
        # Drop the first one, if there is no match on the first position
        if checkfirst and orfs.size>1 and !is_splitter_func.call(orfs.first[0])
          orfs.shift
        end
        orfs
      end

      def get_codon_orfs2 is_splitter_func, is_start_func
        orfs = get_codon_orfs1(is_start_func,false)
        # Check the first one for a start codon
        head = orfs.first
        # Find all start codons and remove leading stop codon
        tail = orfs.find_all { | orf | is_splitter_func.call(orf[1]) }.map { | orf | orf[1..-1] } 
        if is_splitter_func.call(head[0])
          [head] + tail
        else
          tail
        end
      end

      def split codons, func
        list = []
        node = []
        codons.each do | c |
          if func.call(c)
            node.push c
            list.push node if node.size >= @min_size
            node = []
          end
          node.push c
        end
        list
      end

    end

    class OrfEmitter

      # 6-frame ORF emitter for (growing) sequences from the +emit+ 
      # object. Type can be a symbol or a function. Symbols are
      #
      #   :stopstop   All sequences from STOP to STOP codon
      #
      # size control is in nucleotides.
      #
      def initialize emit, type, min_size=30, max_size=nil
        @em = emit
        @type = type
        @min_size = min_size
        @max_size = max_size
      end

      # Concats sequences from the emitter and yields the
      # contained ORFs for every resulting frame (-3..-1, 1..3 )
      def emit_seq
        @em.emit_seq do | part, index, tag, seq |
        end
      end
    end
  end
end
