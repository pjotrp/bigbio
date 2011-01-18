
$: << "../lib"

require 'bigbio'

describe Bio::Big::FastaEmitter, "when using the emitter" do
  include Bio::Big

  it "should emit small parts" do
    s = ""
    FastaEmitter.new("test/data/fasta/nt.fa",10).emit_seq do | part, index, tag, seq |
      # p [index, part, tag, seq]
      s += seq
      if index == 95 and part == :tail
        s.should == "ATCATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGTCTTCAGTGTACAGTATCAAGGGCTCGATCTGCGGTGGATGAGACATCAGATTCAGGAGCTTTTCAAAGAACTGCATCGACATCCGTAACTTCGTTTCAAAAGATTCCAATTCTCAGTTTCAGCTGAATCTGGTAGATACCATCTTTACATATCGTATGCTTGTCATGGGCTTCTAGATGCCTTTCATACTTAAAGATCAAAGGACTTGACGATGCAATAAGCTTCTCGTCTGTAAAACCC"
      end
      s = "" if part == :tail
    end
  end

  it "should emit large parts" do 
    FastaEmitter.new("test/data/fasta/nt.fa").emit_seq do | part, index, tag, seq |
      # p [index, part, tag, seq]
      if index == 95
        seq.should == "ATCATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGTCTTCAGTGTACAGTATCAAGGGCTCGATCTGCGGTGGATGAGACATCAGATTCAGGAGCTTTTCAAAGAACTGCATCGACATCCGTAACTTCGTTTCAAAAGATTCCAATTCTCAGTTTCAGCTGAATCTGGTAGATACCATCTTTACATATCGTATGCTTGTCATGGGCTTCTAGATGCCTTTCATACTTAAAGATCAAAGGACTTGACGATGCAATAAGCTTCTCGTCTGTAAAACCC"
      end
    end
  end
end

describe Bio::Big::FrameState, "when using the FrameState" do

  include Bio::Big

  it "should grow with sequences in frame 1 and return codons" do
    fr = FrameState.new
    fr.seq.should == ''
    fr.add "agtcatc"
    fr.pos.should == 0
    fr.added_codons.should == ['AGT','CAT']
    fr.add "agtcat"
    fr.pos.should == 7
    fr.added_codons.should == ['CAG','TCA']
    fr.add "agt"
    fr.pos.should == 13
    fr.added_codons.should == ['TAG']
  end

  it "should grow with sequences in frame 3 and return codons" do
    fr = FrameState.new
    fr.add "tcatc"
    fr.pos.should == 0
    fr.added_codons.should == ['TCA']
    fr.add "agtcat"
    fr.pos.should == 5
    fr.added_codons.should == ['TCA','GTC']
    fr.add "agt"
    fr.added_codons.should == ['ATA']
  end

  it "should find an ORF in" do
    fr = FrameState.new "atg"
    fr.stopstop?.should == false
    fr = FrameState.new "atggattaaatgtaa"
    fr.stopstop?.should == true
    fr.start.should == 6
    fr.stop.should == 6+6
    fr.fetch.should == "TAAATGTAA"
    fr.fetch.should == nil
  end
  it "should find two ORFs in" do
    fr = FrameState.new "atggattaaatgtaatgttgttaa"
    fr.hasorf?.should == true
    fr.fetch.should == "TAAATGTAA"
    fr.fetch.should == "TAATGTTGTTAA"
    fr.fetch.should == nil
  end
end

describe Bio::Big::OrfEmitter, "when using the ORF emitter" do
  include Bio::Big

  it "should emit STOP-STOP ORFs in all frames" do
    f = FastaEmitter.new("test/data/fasta/nt.fa")
    OrfEmitter.new(f,:stopstop)::emit_seq do | frame, index, tag, seq |
      p [index, tag ] # , seq]
    end
  end
  if false
  it "should emit START-STOP ORFs in all frames"
  it "should emit ORFs on any filter"
  it "should emit ORFs using a minimum size"
  it "should emit ORFs with adjoining sequences"
  end
end
