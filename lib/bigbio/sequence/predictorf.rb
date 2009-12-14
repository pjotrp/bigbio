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
    @seq[@start..@stop-1]
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
  def initialize num, type, id, descr, nt, frame, start, aa
    @id = id +'_'+(num+1).to_s
    # ---- adjust start to match frame
    start += frame.abs-1
    # ---- stop should not go beyond sequence
    stop = start + aa.size * 3
    if stop > nt.size
      stop = nt.size
    end
    # ---- if frame < 0 it should reverse
    nt = nt.reverse if frame < 0
    # p [start, stop, stop-start]
    # p nt
    fr = frame.to_s
    fr = '+'+fr if frame > 0
    @descr = "[#{type} #{fr} #{start} - #{stop}; #{stop-start}/#{nt.size}] " + descr
    @nt = ORFnucleotides.new(nt, start, stop)
    @frame = frame
    @aa = ORFaminoacids.new(aa)
  end

  def <=> orf
    orf.aa.seq.size <=> aa.seq.size
  end

  def to_fastarec
    aa = FastaRecord.new(@id,@descr,@aa.seq)
    nt = FastaRecord.new(@id,@descr,@nt.seq)
    FastaPairedRecord.new(nt,aa)
  end
end

class PredictORF

  def initialize id, descr, seq, trn_table
    @id        = id
    @descr     = descr
    @seq       = seq.gsub(/\s/,'')
    @trn_table = trn_table
  end

  # Return a list of predicted ORFs with :minsize AA's. The ORF's
  # are between STOP codons (so sequences without START are included)
  def stopstop minsize=30
    type = "XX"
    orfs = []
    translate = Nucleotide::Translate.new(@trn_table)
    aa_frames = translate.aa_frames(@seq)
    num = 0
    aa_frames.each do | aa_frame |
      frame = aa_frame[:frame]
      aa = aa_frame[:sequence]
      aa_start = 0
      aa.split(/\*/).each do | candidate |
        if candidate.size >= minsize and candidate.size > 0
          orf = ORF.new(num,type,@id,@descr,@seq,frame,aa_start*3,candidate)
          orfs.push orf
          num += 1
        end
        aa_start += candidate.size + 1
      end
    end
    orfs.sort
  end

  # Return a list of predicted ORFs with :minsize AA's. The ORF's
  # are between START and STOP codons (ATG, TTG, CTG and AUG, UUG and CUG for
  # now, a later version should use the EMBOSS translation table).
  def startstop minsize=30
    stopstop(minsize).find_all { | orf | 
      codon1= orf.nt.seq[0..2].upcase
      ['ATG','TTG','CTG','AUG','UUG','CUG'].index(codon1) != nil
    }
  end 

  # Return the longest ORF that has a START codon (see +startstop+)
  # Returns nil if none is found
  def longest_startstop
    startstop(0).first
  end
    
end
