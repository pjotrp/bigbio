
require 'set'

module Bio
  module Big

    module FrameCodonHelpers

      STOP_CODONS = Set.new(%w{TAG TAA TGA UAG UAA UGA})
      START_CODONS = Set.new(%w{ATG AUG})

      # Track sequence position in parent sequence
      module TrackSequenceTrait
        attr_accessor :track_sequence_pos 
        def TrackSequenceTrait.update_sequence_pos orfs, seq_pos
          orfs.each { | orf | orf.track_sequence_pos = seq_pos + orf.pos }
        end
      end

      # Functions that move a frame forward, or backward, 
      # creating new short frames.
      module MoveShortFrame
        def shortframe_right seq
          fr = self
          p fr
        end
      end

      class FrameCodonSequence 
        include Enumerable
        include TrackSequenceTrait
        attr_reader :pos     # codon position in short parent sequence
        attr_reader :codons
        def initialize seq, pos=0
          if seq.kind_of?(String)
            @codons = seq.upcase.scan(/(\w\w\w)/).flatten
          else
            @codons = seq
          end
          @pos = pos
        end
        def size
          @codons.size
        end
        def rpos
          pos + size
        end
        def [] index
          @codons[index]
        end
        def shift
          list = @codons
          list.shift
          FrameCodonSequence.new(list,@pos+1)
        end
        def to_seq
          @codons.join
        end
        def each
          @codons.each { | c| yield c }
        end
      end
    end # FrameCodonHelpers

    # The short frame uses the simplest concept to find ORFs. The sequence is
    # immutable, always forward and in frame 0. That makes it easy to reason.
    # It also return all ORF's in one go, with the left/right locations.

    class ShortFrameState
      include FrameCodonHelpers
      include MoveShortFrame

      def initialize seq, min_size = 30
        # @seq = seq.upcase  
        @min_size_codons = (min_size/3).to_i
        @codons = FrameCodonSequence.new(seq)
      end

      # Return a list of ORFs delimited by STOP codons. 
      def get_stopstop_orfs 
        get_codon_orfs1(Proc.new { | codon | STOP_CODONS.include?(codon) },false,true)
      end

      # Return a list of ORFs delimited by START-STOP codons
      def get_startstop_orfs 
        get_codon_orfs2(
                 Proc.new { | codon | STOP_CODONS.include?(codon) },
                 Proc.new { | codon | START_CODONS.include?(codon) })
      end

      # Splitter for one delimiter function. +include_leftmost+ decides
      # the first sequence is returned when incomplete. +strip_leading+
      # is used to remove the shared codon with the last sequence.
      #
      def get_codon_orfs1 splitter_func,do_include_leftmost,do_strip_leading 
        orfs = split(@codons,splitter_func)
        # Drop the first sequence, if there is no match on the first position
        orfs.shift if !do_include_leftmost and orfs.size>1 and !splitter_func.call(orfs.first[0])
        orfs.map { |codons| 
          codons = codons.shift if do_strip_leading and splitter_func.call(codons[0])
          codons
        }
      end

      # Splitter for two delimeter functions
      def get_codon_orfs2 splitter_func, start_func
        orfs = get_codon_orfs1(splitter_func,true,true)
        orfs.find_all { | orf | start_func.call(orf[0]) }
      end

      # Return list of codon sequences, split on the +is_splitter+ 
      # function.
      #
      def split codons, is_splitter_func
        list = []
        node = []
        codons.each_with_index do | c, pos |
          # p [c,pos]
          if is_splitter_func.call(c)
            node.push c
            size = node.size
            list.push FrameCodonSequence.new(node,pos+1-size) if size > @min_size_codons
            node = []
          end
          node.push c
        end
        list
      end

    end

    module GetFrame

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
          # p [part, seq]
          break
        end
      end
    end
  end
end
