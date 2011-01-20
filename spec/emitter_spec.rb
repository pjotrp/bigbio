
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

describe Bio::Big::ShortFrameState, "when using the ShortFrameState" do

  include Bio::Big
  
  it "should find an ORF" do
    fr = ShortFrameState.new "atggattaaatgtaatggatttaatgtaaa",0
    orfs = fr.get_stopstop_orfs
    orfs.map{ | orf | orf.pos }.should == [ 3, 5 ]
    orfs.map{ | orf | orf.to_seq }.should == ["ATGTAA", "TGGATTTAA"]
    orfs = fr.get_startstop_orfs
    orfs.map{ | orf | orf.pos }.should == [ 0, 3 ]
    orfs.map{ | orf | orf.to_seq }.should == ["ATGGATTAA","ATGTAA"]
  end
  it "should handle min_size" do
    fr = ShortFrameState.new "atggattaaatgtaatggatttaatgtaaa",9
    orfs = fr.get_stopstop_orfs
    orfs.map{ | orf | orf.pos }.should == [ 5 ]
    orfs.map{ | orf | orf.to_seq }.should == [ "TGGATTTAA"]
    fr.get_startstop_orfs.should == []
  end
  it "should find ORFs in" do
    fr = ShortFrameState.new "atgttttaaatgtaatgttgttaaatgttttaaatgtaatgttgttaa",0
    orfs = fr.get_stopstop_orfs
    orfs.map{ | orf | orf.to_seq }.should == ["ATGTTTTAA", "ATGTAA", "ATGTTTTAA", "ATGTAA"]
    orfs.map{ | orf | orf.pos }.should == [ 0, 3, 8, 11]
    orfs = fr.get_startstop_orfs
    orfs.map{ | orf | orf.to_seq }.should == ["ATGTTTTAA", "ATGTAA", "ATGTTTTAA", "ATGTAA"]
    orfs.map{ | orf | orf.pos }.should == [ 0, 3, 8, 11]
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
