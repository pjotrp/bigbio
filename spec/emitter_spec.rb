
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
    fr.get_stopstop_orfs.should == ["ATGTAA", "TGGATTTAA"]
    fr.get_startstop_orfs.should == ["ATGGATTAA", "ATGTAA"]
  end
  it "should handle min_size" do
    fr = ShortFrameState.new "atggattaaatgtaatggatttaatgtaaa",9
    fr.get_stopstop_orfs.should == ["TGGATTTAA"]
    fr.get_startstop_orfs.should == []
  end


end

describe Bio::Big::FrameState, "when using the FrameState" do

  include Bio::Big

  if false
  it "should grow with sequences in frame 1 and return codons" do
    fr = FrameState.new
    fr.seq.should == ''
    fr.append "agtcatc"
    fr.pos.should == 0
    fr.added_codons.should == ['AGT','CAT']
    fr.append "agtcat"
    fr.pos.should == 7
    fr.added_codons.should == ['CAG','TCA']
    fr.append "agtt"
    fr.pos.should == 13
    fr.added_codons.should == ['TAG']
    fr.append "agt"
    fr.added_codons.should == ['TTA']
  end

  it "should grow with sequences in frame 3 and return codons" do
    fr = FrameState.new
    fr.append "tcatc"
    fr.pos.should == 0
    fr.added_codons.should == ['TCA']
    fr.append "agtcat"
    fr.pos.should == 5
    fr.added_codons.should == ['TCA','GTC']
    fr.append "agt"
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
  it "should find two STOP_STOP ORFs in" do
    fr = FrameState.new "atggattaaatgtaatgttgttaa"
    fr.hasorf?.should == true
    fr.fetch.should == "TAAATGTAA"
    fr.fetch.should == "TAATGTTGTTAA"
    fr.fetch.should == nil
  end
  it "should find two ORFs in" do
    fr = FrameState.new "atgttttaaatgtaatgttgttaa", :startstop
    fr.hasorf?.should == true
    fr.fetch.should == "ATGTTTTAA"
    fr.fetch.should == "ATGTAA"
    fr.fetch.should == nil
  end
  it "should find four ORFs in" do
    fr = FrameState.new "atgttttaaatgtaatgttgttaa", :startstop
    fr.hasorf?.should == true
    fr.fetch.should == "ATGTTTTAA"
    fr.fetch.should == "ATGTAA"
    fr.append "atgttttaaatgtaatgttgttaa"
    fr.hasorf?.should == true
    fr.fetch.should == "ATGTTTTAA"
    fr.fetch.should == "ATGTAA"
    fr.fetch.should == nil
  end
  end
end

describe Bio::Big::ReversedFrameState, "when using the ReversedFrameState" do

  if false
  include Bio::Big

  it "should grow with sequences in frame 1 and return codons" do
    fr = ReversedFrameState.new
    fr.seq.should == ''
    fr.prepend "tga"
    # 5' tga 3'
    fr.pos.should == 3
    fr.added_codons.should == ['TGA']
    fr.prepend "tactga"
    # 5' tactga tga 3'
    fr.pos.should == 6
    fr.added_codons.should == ['TGA','TAC']
    fr.prepend "ctactga"
    # 5' c|tactga tactga tga 3'
    fr.pos.should == 7
    fr.c_pos.should == 7
    fr.added_codons.should == ['TGA','TAC']
    fr.prepend "tcagtc"
    # 5' tcagtc c|tactga tactga tga 3'
    fr.pos.should == 6
    fr.c_pos.should == 7
    fr.added_codons.should == ['TCC','CAG']
  end

  it "should find an ORF in" do
    fr = ReversedFrameState.new "atgatg"
    fr.stopstop?.should == false
    # 5' tagtttaaatag 3'
    fr = ReversedFrameState.new "atggattaaatgtaa"
    fr.added_codons[0].should == "TAA"
    fr.stopstop?.should == true
    fr.fetch.should == "TAAATGTAA"
    fr.fetch.should == nil
  end
  it "should find two ORFs in" do
    fr = ReversedFrameState.new "atggattaaatgtaatgttgttaa"
    fr.hasorf?.should == true
    fr.fetch.should == "TAATGTTGTTAA"
    p [fr,fr.added_codons]
    fr.fetch.should == "TAAATGTAA"
    fr.fetch.should == nil
  end
  it "should find two ORFs in" do
    fr = ReversedFrameState.new "atgttttaaatgtaatgttgttaa", :startstop
    fr.hasorf?.should == true
    fr.fetch.should == "ATGTAATGTTGTTAA"
    fr.fetch.should == "ATGTTTTAA"
    fr.fetch.should == nil
  end
  it "should find four ORFs in" do
    fr = ReversedFrameState.new "atgttttaaatgtaatgttgttaa", :startstop
    fr.hasorf?.should == true
    # p fr
    # p fr.added_codons
    fr.fetch.should == "ATGTAATGTTGTTAA"
    fr.fetch.should == "ATGTTTTAA"
    fr.hasorf?.should == false
    # p fr
    fr.prepend "atgttttaaatgtaatgttgttaa"
    fr.hasorf?.should == true
    # p fr
    # p fr.added_codons
    fr.fetch.should == "ATGTAATGTTGTTAA"
    fr.fetch.should == "ATGTTTTAA"
    fr.fetch.should == nil
  end
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
