
require 'set'

module Bio
  module Big

    STOP_CODONS = Set.new(%w{TAG TAA TGA UAG UAA UGA})
    START_CODONS = Set.new(%w{ATG AUG})

    module FrameStateHelpers
      # Finds the next type location, starting from pos-retrack. 
      #
      # True if found, otherwise false
      def hasorf?
        case @type
          when :stopstop then stopstop?
          when :startstop then startstop?
          else
            raise "Unknown type to act on #{@type}"
        end
      end

      def stopstop?
        found?(Proc.new { | codon | STOP_CODONS.include?(codon) }, 
               Proc.new { | codon | STOP_CODONS.include?(codon) })
      end

      def startstop?
        found?(Proc.new { | codon | START_CODONS.include?(codon) }, 
               Proc.new { | codon | STOP_CODONS.include?(codon) })
      end
    end

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

      def get_codon_orfs1 func, checkfirst=true
        orfs = split(@codons,func)
        # Drop the first one, if there is no match on the first position
        if checkfirst and orfs.size>1 and !func.call(orfs.first[0])
          orfs.shift
        end
        orfs
      end

      def get_codon_orfs2 func1, func2
        orfs = get_codon_orfs1(func2,false)
        p orfs
        # Check the first one for a start codon
        head = orfs.first
        tail = orfs.find_all { | orf | func1.call(orf[1]) }.map { | orf | orf[1..-1] } 
        if func1.call(head[0])
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

      def found? func1, func2
        codons = added_codons
        codon1 = 0
        if @start == nil
          # look for first STOP codon
          codons.each_with_index { | codon, idx | 
            if func1.call(codon)
              codon1 = idx
              @start = idx * 3 + @c_pos
              break
            end
          }
        end
        if @start != nil and @stop == nil
          # look for 2nd STOP codon
          codons[codon1+1..-1].each_with_index { | codon, idx |
            if func2.call(codon)
              # p [idx,codon]
              @stop = (codon1 + 1 + idx)*3 + @c_pos 
              break
            end
          }
        end
        return (@start!=nil and @stop!=nil)
      end

    end

    class FrameState

      include FrameStateHelpers

      attr_reader :seq, :pos, :type, :start, :stop

      # Keeps track of a frame by adding partial sequences and scanning
      # for ORFs. The sequence should be in frame to keep reasoning 
      # easy. If a find is made the sequence should be reset.
      # frame (0..2)
      def initialize seq = '', type=:stopstop 
        @type = type
        @seq = seq.upcase
        @pos = 0
        @c_pos = 0
        @start = nil  # keep track of first find
        @stop  = nil 
      end

      def append seq
        if hasorf?
          raise "You need to fetch ORFs before appending!"
        end
        @pos = @seq.size 
        @c_pos = @pos - @pos % 3 # round to nearest CODON
        @seq += seq.upcase
      end

      # Fetch the ORF and reset state to start at new part of sequence
      def fetch
        if hasorf?
          orf = @seq[@start..@stop+2]
          @seq = @seq[@stop..-1]  # Retain last codon!
          @start = nil
          @stop = nil
          @pos = 0
          @c_pos = 0
          orf
        else
          nil
        end
      end

      def found? func1, func2
        codons = added_codons
        codon1 = 0
        if @start == nil
          # look for first STOP codon
          codons.each_with_index { | codon, idx | 
            if func1.call(codon)
              codon1 = idx
              @start = idx * 3 + @c_pos
              break
            end
          }
        end
        if @start != nil and @stop == nil
          # look for 2nd STOP codon
          codons[codon1+1..-1].each_with_index { | codon, idx |
            if func2.call(codon)
              # p [idx,codon]
              @stop = (codon1 + 1 + idx)*3 + @c_pos 
              break
            end
          }
        end
        return (@start!=nil and @stop!=nil)
      end

      # Return a list of codons to check
      def added_codons 
        seq1 = @seq[@c_pos..-1]
        seq2 = seq1.scan(/(\w\w\w)/).flatten
        seq2
      end

    end

    class ReversedFrameState
      include FrameStateHelpers

      attr_reader :seq, :pos, :c_pos, :type, :start, :stop

      # Very similar to FrameState, but prepend sequences, rather that
      # adding at the end. A seperate algorithm is needed to keep 
      # speed up. In this edition we look for the STOP signal first.
      #
      # Keeps track of a frame by adding partial sequences and scanning for
      # ORFs. The sequence should be in (right-side) frame to keep reasoning
      # easy. If a find is made, the sequence should be reset. frame (0..2)
      def initialize seq = '', type=:stopstop 
        @type = type
        @seq = seq.upcase
        @pos = seq.size
        @c_pos = @pos + @seq.size % 3 # round to nearest CODON
        @start = nil  # keep track of first find
        @stop  = nil  # keep track of last find
      end

      def prepend seq
        if hasorf?
          raise "You need to fetch ORFs before prepending!"
        end
        @pos = seq.size
        @c_pos = @pos + @seq.size % 3 # round to nearest CODON
        @seq = seq.upcase + @seq
      end

      # Fetch the ORF and reset state to start at new part of sequence
      def fetch
        if hasorf?
          len = @seq.size
          p1 = len - @start - 3  # start counts from right
          p2 = len - @stop - 1   # stop counts from right
          # ---------------p1-----------p2--------cpos--- (len)
          #               orf----------orf
          # seq--------------seq
          # --------------pos/c_pos
          orf = @seq[p1..p2]
          @seq = @seq[0..p1+2]
          @pos = p1
          @c_pos = @pos + @seq.size % 3 # round to nearest CODON
          @start = nil
          @stop = nil
          orf
        else
          nil
        end
      end

      def found? func1, func2
        codons = added_codons
        codon1 = 0
        if @stop == nil
          # look for right/stop codon
          codons.each_with_index { | codon, idx | 
            if func2.call(codon)
              codon1 = idx
              @stop = idx * 3  # from the end!
              break
            end
          }
        end
        if @stop != nil and @start == nil
          # look for left/start codon
          codons[codon1+1..-1].each_with_index { | codon, idx |
            if func1.call(codon)
              # p [idx,codon]
              @start = (codon1+1 + idx)*3 
              break
            end
          }
        end
        return (@start!=nil and @stop!=nil)
      end

      # Return a list of codons to check
      def added_codons 
        seq1 = @seq[@c_pos%3..@c_pos]
        seq2 = seq1.scan(/(\w\w\w)/).flatten.reverse
        seq2
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
