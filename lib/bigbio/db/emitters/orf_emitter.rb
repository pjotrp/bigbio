
module Bio
  module Big
    class OrfEmitter

      # 6-frame ORF emitter for (growing) sequences from the +emit+ 
      # object. Type can be a symbol or a funtion (NYI). Symbols are
      #
      #   :stopstop   All sequences from STOP to STOP codon
      def initialize emit, type, min_size, max_size
        @em = emit
        @type = type
        @min_size = min_size
        @max_size = max_size
      end

      def emit_seq
      end
    end
  end
end
