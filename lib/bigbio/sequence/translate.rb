
require 'biolib/emboss'

module Nucleotide

  class Translate

      # Table can be either an id (integer) or a Biolib::Emboss TrnTable

      def initialize table
        if table.kind_of? Numeric
          @trn_table = Biolib::Emboss.ajTrnNewI(table)
        else
          @trn_table = table
        end
      end

      # Return all six reading frames as an Array - ordered as
      # frames [1,2,3,-1,-2,-3] with as tuples [frame, AAsequence]

      def aa_frames seq
        res = []
        # remove white space
        seq = seq.gsub(/\s/,'')
        ajpseq     = Biolib::Emboss.ajSeqNewNameC(seq,"Test sequence")
        [1,2,3,-1,-2,-3].each do | frame |
          ajpseqt  = Biolib::Emboss.ajTrnSeqOrig(@trn_table,ajpseq,frame)
          aa       = Biolib::Emboss.ajSeqGetSeqCopyC(ajpseqt)
          res.push({:frame => frame, :sequence => aa})
        end
        res
      end
    end

end

