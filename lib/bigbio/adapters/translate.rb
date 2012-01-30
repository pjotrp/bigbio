# TranslationAdapter will translate using EMBOSS, or BioRuby
# when the first is not available

module Bio
  module Big
    module TranslationAdapter

      VALID_FRAME_VALUES = [ 0, -1, -2, -3, 1, 2, 3 ]
      def self.translation_table num
        if Environment.instance.biolib
          Biolib::Emboss.ajTrnNewI(num)
        end
      end

      # Precompile sequence for EMBOSS
      def self.pre_translate seq,label
        if Environment.instance.biolib
          Biolib::Emboss.ajSeqNewNameC(seq,"Test sequence")
        else
          nil
        end
      end

      # Translate using frame (pre_seq is only used for EMBOSS)
      # 
      # Valid frame values are 0,1,2,3 and -1,-2,-3, where 0 and 1 are the
      # standard reading frame. The negative values translate the reverse
      # complement of the strand.
      def self.translate trn_table, frame, seq, pre_seq = nil
        raise "Illegal frame #{frame}" if VALID_FRAME_VALUES.index(frame) == nil
        frame = 1 if frame == 0
        if Environment.instance.biolib
          # Using EMBOSS for translation
          ajpseq = pre_seq
          if not pre_seq
            ajpseq = Biolib::Emboss.ajSeqNewNameC(seq,"Test sequence")
          end
          ajpseqt  = Biolib::Emboss.ajTrnSeqOrig(trn_table,ajpseq,frame)
          Biolib::Emboss.ajSeqGetSeqCopyC(ajpseqt)
        else
          # Using BioRuby for translation
          ntseq = if frame > 0 
            Bio::Sequence::NA.new(seq[frame-1..-1])
          else
            # This to match EMBOSS frames
            rframe =
              case frame
                when -2 
                  -3 
                when -3
                  -2
                else 
                  -1
              end
            Bio::Sequence::NA.new(seq[0..rframe]).reverse_complement
          end
          # pp ntseq
          ntseq.translate.to_s
        end
      end
    end
  end
end
