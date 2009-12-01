# ORF predictor class
#

require 'bigbio/sequence/translate'

class ORFsequence
  attr_accessor :seq
  def initialize sequence
    @seq = sequence
  end
end

# Helper class for storing ORF information
class ORFnucleotides < ORFsequence
  attr_reader :start, :stop
  def initialize sequence, start, stop
    super(sequence)
    @start = start
    @stop = stop
  end

  def seq
    @seq[@start..@stop]
  end

  def fullseq
    @seq
  end
end

# Helper class for storing ORF information
class ORFaminoacids < ORFsequence
end


class ORF
  attr_reader :id, :descr, :nt, :aa, :frame
  def initialize id, descr, nt, frame, start, aa
    @id = id
    @descr = descr
    @nt = ORFnucleotides.new(nt, start, start + aa.size * 3)
    @frame = frame
    @aa = ORFaminoacids.new(aa)
  end

  def <=> orf
    orf.aa.seq.size <=> aa.seq.size
  end

end

class PredictORF

  def initialize id, descr, seq, trn_table
    @id        = id
    @descr     = descr
    @seq       = seq.gsub(/\s/,'')
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
      aa_start = 0
      aa.split(/\*/).each do | candidate |
        if candidate.size >= minsize
          orf = ORF.new(@id,@descr,@seq,frame,aa_start*3,candidate)
          orfs.push orf
        end
        aa_start += candidate.size + 1
      end
    end
    orfs.sort
  end

end
