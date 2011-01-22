
require 'set'

module Bio
  module Big

    module FrameCodonHelpers

      STOP_CODONS = Set.new(%w{TAG TAA TGA UAG UAA UGA})
      START_CODONS = Set.new(%w{ATG AUG})

      # Track sequence position in parent sequence (in nucleotides)
      module TrackSequenceTrait
        attr_accessor :track_ntseq_pos 
        def TrackSequenceTrait.update_sequence_pos orfs, ntseq_pos
          orfs.each { | orf | orf.track_ntseq_pos = ntseq_pos + orf.pos*3 }
          orfs
        end
        def TrackSequenceTrait.update_reversed_sequence_pos orfs, ntseq_pos
          # is the same
          orfs.each { | orf | orf.track_ntseq_pos = ntseq_pos + orf.pos*3 }
          orfs
        end
      end

      # Functions that move a frame forward, or backward, 
      # creating new short frames.
      module CreateShortFrame

        def CreateShortFrame.create_right fr,orfs,rseq
          seq = fr.seq
          ntseq_pos = fr.ntseq_pos
          remove = if orfs.size > 0
            orfs.last.rpos*3
          else 
            0
          end
          ntseq_pos += remove
          nseq = seq[remove..-1] + rseq
          ShortFrameState.new nseq,ntseq_pos,fr.min_size_codons*3
        end

        def CreateShortFrame.create_left fr,orfs,nseq
          # Reversed (real locations on contig):
          #
          # |  3                21  B |
          # ttaaatgtaatttaggtaaatttat atgtaaattaggta (reversed)
          # ...^--============xxx^=======xxx
          #       ^                     ^
          # Actual feed:
          #
          # s2=              s1=
          # "atggattaaatgta" "tatttaaatggatttaatgtaaatt"
          #  ......xxx=====   ~===xx^============--^...                               
          #  0  1  2  3        0  1  2  3 
          seq1 = fr.seq             # original sequence
          len1 = seq1.size
          ntseq_pos1 = fr.ntseq_pos # right side of seq (|)
          bridge = len1 % 3    # chomp left side (B)
          remove = if orfs.size > 0
            len1 - bridge - (orfs.first.pos)*3 + 1
          else 
            0
          end
          ntseq_pos2 = ntseq_pos1+remove-1  # pos against main contig
          seq2 = nseq + seq1[0..(len1-remove)]
          ShortReversedFrameState.new seq2,ntseq_pos2,fr.min_size_codons*3
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
      attr_reader :seq, :ntseq_pos, :min_size_codons, :codons

      def initialize seq, ntseq_pos, ntmin_size
        # @seq = seq.upcase  
        @seq = seq
        @min_size_codons = if ntmin_size > 3
                             (ntmin_size/3).to_i
                           else
                             2  # otherwise we get single STOP codons
                           end
       
        @codons = FrameCodonSequence.new(seq,ntseq_pos)
        @ntseq_pos = ntseq_pos # nucleotides
        # @codons.track_sequence_pos = seq_pos
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
      def get_codon_orfs1 splitter_func,do_include_leftmost_orf,do_strip_leading_codon
        orfs = split(@codons,splitter_func)
        return [] if orfs.size == 0
        # Drop the first sequence, if there is no match on the first position
        orfs.shift if !do_include_leftmost_orf and !splitter_func.call(orfs.first[0])
        orfs = orfs.map { |codons| 
          codons = codons.shift if do_strip_leading_codon and splitter_func.call(codons[0])
          codons
        }
        if @reversed == nil
          TrackSequenceTrait.update_sequence_pos(orfs,@ntseq_pos) # nail against parent
        else
          TrackSequenceTrait.update_reversed_sequence_pos(orfs,@ntseq_pos) # nail against parent
        end
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
            # p node
            list.push FrameCodonSequence.new(node,pos+1-size) if size > @min_size_codons
            node = []
          end
          node.push c  # always push boundary codon
        end
        list
      end

    end

    # This is the reversed version, which is rather the same as the forward,
    # though the tracked ntseq_pos should be seen from the end of the sequence,
    # as we are emmiting sequences from the end(!) Also we need to make sure
    # the sequence is always in frame (from the left).
    class ShortReversedFrameState < ShortFrameState
      attr_accessor :reversed
  
      def initialize seq, ntseq_pos, ntmin_size
        @reversed = true
        chop = seq.size % 3 # align on codons
        super seq[chop..-1],ntseq_pos,ntmin_size
        @seq = seq # but record full seq
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
          # p [part, seq]
          break
        end
      end
    end
  end
end
