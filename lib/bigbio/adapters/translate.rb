# TranslationAdapter will translate using EMBOSS, or BioRuby
# when the first is not available

module Bio
  module Big
    module TranslationAdapter
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

      # Translate
      def self.translate trn_table, frame, seq, pre_seq = nil
        if Environment.instance.biolib
          ajpseq = pre_seq
          if not pre_seq
            ajpseq = Biolib::Emboss.ajSeqNewNameC(seq,"Test sequence")
          end
          ajpseqt  = Biolib::Emboss.ajTrnSeqOrig(trn_table,ajpseq,frame)
          Biolib::Emboss.ajSeqGetSeqCopyC(ajpseqt)
        else
          ntseq = if frame > 0 
            Bio::Sequence::NA.new(seq[frame-1..-1])
          else
            Bio::Sequence::NA.new(seq.reverse[-frame-1..-1])
          end
          # ntseq
          ntseq.translate.to_s
        end
      end
    end
  end
end
