

module Nucleotide

  module TranslationTable
  end

  class Translate

    include Bio::Big::TranslationAdapter

      # Table can be either an id (integer) or a Biolib::Emboss TrnTable

      def initialize table
        if table.kind_of? Numeric
          @trn_table = translation_table(table)
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
        [1,2,3,-1,-2,-3].each do | frame |
          aa = translate(@trn_table,frame,seq)
          res.push({:frame => frame, :sequence => aa})
        end
        res
      end
    end

end

