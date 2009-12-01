# ORF predictor class
#

require 'bigbio/sequence/translate'

# Helper class for storing ORF information
class ORFnucleotides
end

# Helper class for storing ORF information
class ORFaminoacids
end


class ORF
  attr_reader :nt_sequence, :aa_sequence, :frame
  def initialize nt, frame, aa
    @nt = nt
    @frame = frame
    @aa = aa
  end
end

class PredictORF

  def initialize id, descr, seq, trn_table
    @id        = id
    @descr     = descr
    @seq       = seq
    @trn_table = trn_table
  end

  # Return a list of predicted ORFs with :minsize AA's
  def stopstop minsize=30
    orfs = []
    translate = Nucleotide::Translate.new(@trn_table)
    aa_frames    = translate.aa_frames(@seq)
    aa_frames.each do | aa_frame |
      frame = aa_frame[:frame]
      aa = aa_frame[:sequence]
      aa.split(/\*/).each do | candidate |
        p candidate
        if candidate.size >= minsize
          orf = ORF.new(@seq,frame,candidate)
          orfs.push orf
        end
      end
    end
    orfs
  end

end
