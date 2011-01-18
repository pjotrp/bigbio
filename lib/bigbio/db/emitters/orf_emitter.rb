
require 'set'

module Bio
  module Big

    STOP_CODONS = Set.new(%w{TAG TAA TGA UAG UAA UGA})
    START_CODONS = Set.new(%w{ATG AUG})

    class FrameState

      attr_reader :seq, :pos, :type, :start, :stop

      # Keeps track of a frame by adding partial sequences and scanning
      # for ORFs. The sequence should be in frame to keep reasoning 
      # easy. If a find is make the sequence should be reset.
      # frame (0..2)
      def initialize seq = '', type=:stopstop 
        @type = type
        @seq = seq.upcase
        @pos = 0
        @c_pos = 0
        @start = nil  # keep track of first find
        @stop  = nil 
      end

      def add seq
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
