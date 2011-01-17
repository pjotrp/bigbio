

module Bio
  module Big
    module TranslationAdapter
      def translation_table num
        if Environment.instance.biolib
          Biolib::Emboss.ajTrnNewI(num)
        end
      end

      def translate trn_table, frame, seq
        if Environment.instance.biolib
          ajpseq   = Biolib::Emboss.ajSeqNewNameC(seq,"Test sequence")
          ajpseqt  = Biolib::Emboss.ajTrnSeqOrig(trn_table,ajpseq,frame)
          Biolib::Emboss.ajSeqGetSeqCopyC(ajpseqt)
        else
          ntseq = if frame > 0 
            Bio::Sequence::NA.new(seq[frame-1..-1])
          else
            Bio::Sequence::NA.new(seq.reverse[-frame-1..-1])
          end
          ntseq.translate.to_s
        end
      end
    end
  end
end
