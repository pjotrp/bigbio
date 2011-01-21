
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
    fr = ShortFrameState.new "atggattaaatgtaatggatttaatgtaaa",0,0
    orfs = fr.get_stopstop_orfs
    orfs.map{ | orf | orf.pos }.should == [ 3, 5 ]
    orfs.map{ | orf | orf.to_seq }.should == ["ATGTAA", "TGGATTTAA"]
    orfs = fr.get_startstop_orfs
    orfs.map{ | orf | orf.pos }.should == [ 0, 3 ]
    orfs.map{ | orf | orf.to_seq }.should == ["ATGGATTAA","ATGTAA"]
  end
  it "should handle min_size" do
    fr = ShortFrameState.new "atggattaaatgtaatggatttaatgtaaa",0,9
    orfs = fr.get_stopstop_orfs
    orfs.map{ | orf | orf.to_seq }.should == [ "TGGATTTAA"]
    orfs.map{ | orf | orf.pos }.should == [ 5 ]
    fr.get_startstop_orfs.should == []
  end
  it "should find ORFs in" do
    fr = ShortFrameState.new "atgttttaaatgtaatgttgttaaatgttttaaatgtaatgttgttaa",0,0
    orfs = fr.get_stopstop_orfs
    orfs.map{ | orf | orf.to_seq }.should ==  ["ATGTAA", "TGTTGTTAA", "ATGTTTTAA", "ATGTAA", "TGTTGTTAA"]
    orfs.map{ | orf | orf.pos }.should == [3, 5, 8, 11, 13]
    orfs = fr.get_startstop_orfs
    orfs.map{ | orf | orf.to_seq }.should == ["ATGTTTTAA", "ATGTAA", "ATGTTTTAA", "ATGTAA"]
    orfs.map{ | orf | orf.pos }.should == [ 0, 3, 8, 11]
  end

  it "should match results of EMBOSS getorf" do
    s = "AG GTTCGNACGGTCATCGNATNAAGTCTTGNATATCG TAA TTNCGCGTGCCGCCTTCTTTCTCCTTTTTCTCTTTTACTTCTTCATCATCATCTTCTTCTTCTTCCTCTTCGATATTCGTCAGTGTGTGTATTTTG GGG AAAACTTTG TGA GCAAAGAGCGAGAAAATGAGCGGANCGG TAA GAAAATCGCGGATGTGGCTTTCAAAGCTTCAAGGACTATCGATTGGGATGGTATGGC TAA GGTCCTTGTCACAGATGAGGCTCGTAGAG".gsub(/ /,'')
    # >_3 [3 - 167] #0
    # 1st  GTTCGNACGGTCATCGNATNAAGTCTTGNATATCGTAATTNCGCGTGCCGCCTTCTTTCT
    # CCTTTTTCTCTTTTACTTCTTCATCATCATCTTCTTCTTCTTCCTCTTCGATATTCGTCA
    # GTGTGTGTATTTTGGGGAAAACTTTGTGAGCAAAGAGCGAGAAAA
    # >_4 [171 - 179] #0
    # OK   GCGGANCGG
    # >_5 [183 - 239] #0
    # OK   GAAAATCGCGGATGTGGCTTTCAAAGCTTCAAGGACTATCGATTGGGATGGTATGGC
    # >_6 [243 - 257] #0
    # OK-  GGTCCTTGTCACAGA
    # >_7 [261 - 266] #0
    # OK-  GGCTCGTAG
    # >_8 [1 - 270] # 1
    # whole!  AGGTTCGNACGGTCATCGNATNAAGTCTTGNATATCGTAATTNCGCGTGCCGCCTTCTTT
    # CTCCTTTTTCTCTTTTACTTCTTCATCATCATCTTCTTCTTCTTCCTCTTCGATATTCGT
    # CAGTGTGTGTATTTTGGGGAAAACTTTGTGAGCAAAGAGCGAGAAAATGAGCGGANCGGT
    # AAGAAAATCGCGGATGTGGCTTTCAAAGCTTCAAGGACTATCGATTGGGATGGTATGGCT
    # AAGGTCCTTGTCACAGATGAGGCTCGTAGA
    # >_1 [2 - 37] #2
    # 1st-   GGTTCGNACGGTCATCGNATNAAGTCTTGNATATCG
    # >_2 [41 - 148] #2
    # OK-   TTNCGCGTGCCGCCTTCTTTCTCCTTTTTCTCTTTTACTTCTTCATCATCATCTTCTTCT
    # TCTTCCTCTTCGATATTCGTCAGTGTGTGTATTTTGGGGAAAACTTTG
    # >_9 [152 - 271] #2
    # last- CAAAGAGCGAGAAAATGAGCGGANCGGTAAGAAAATCGCGGATGTGGCTTTCAAAGCTT
    # CAAGGACTATCGATTGGGATGGTATGGCTAAGGTCCTTGTCACAGATGAGGCTCGTAGAG

    # Frame 0
    minsize = 0
    fr = ShortFrameState.new s,0,minsize
    orfs = fr.get_stopstop_orfs
    orfs.map{ | orf | orf.to_seq }.should == []

    # Frame 1
    fr = ShortFrameState.new s[1..-1],0,minsize
    os = fr.get_stopstop_orfs
    os.map{ | orf | orf.to_seq }.should == ["TTNCGCGTGCCGCCTTCTTTCTCCTTTTTCTCTTTTACTTCTTCATCATCATCTTCTTCTTCTTCCTCTTCGATATTCGTCAGTGTGTGTATTTTGGGGAAAACTTTGTGA"]
    orfs += os
    # Frame 2
    fr = ShortFrameState.new s[2..-1],0,minsize
    os = fr.get_stopstop_orfs
    os.map{ | orf | orf.to_seq }.should == ["GCGGANCGGTAA", "GAAAATCGCGGATGTGGCTTTCAAAGCTTCAAGGACTATCGATTGGGATGGTATGGCTAA", "GGTCCTTGTCACAGATGA", "GGCTCGTAG"]
    orfs += os
    orfs.size.should == 5

    # >_1 [235 - 270] 
    # Last: ATGGCTAAGGTCCTTGTCACAGATGAGGCTCGTAGA
    # >_2 [167 - 271] 
    # Last: ATGAGCGGANCGGTAAGAAAATCGCGGATGTGGCTTTCAAAGCTTCAAGGACTATCGATT
    # GGGATGGTATGGCTAAGGTCCTTGTCACAGATGAGGCTCGTAGAG

    # Frame 0
    minsize = 0
    fr = ShortFrameState.new s,0,minsize
    orfs = fr.get_startstop_orfs

    # Frame 1
    fr = ShortFrameState.new s[1..-1],0,minsize
    orfs += fr.get_startstop_orfs
    # Frame 2
    fr = ShortFrameState.new s[2..-1],0,minsize
    orfs += fr.get_startstop_orfs
    orfs.map{ | orf | orf.to_seq }.should == []
    orfs.size.should == 0
  end
end

describe Bio::Big::ShortFrameState, "when combining frames" do
  include Bio::Big
  it "should combine a forward frame" do
    s1 = "atggattaaatgtaata"
    s2 = "atggatttaatgtaaa"
    fr = ShortFrameState.new s1,0,0
    fr.ntseq_pos.should == 0
    orfs = fr.get_stopstop_orfs
    orfs.size == 1 # in codons
    fr3 = FrameCodonHelpers::CreateShortFrame.create_right(fr,orfs,s2)
    fr3.ntseq_pos.should == 15
    fr3.codons.to_seq.should == "TAATGGATTTAATGTAAA"
    norfs = fr3.get_stopstop_orfs
    orfs += norfs
    orfs.map{ | orf | orf.to_seq }.should == ["ATGTAA", "TGGATTTAA"]
    orfs.map{ | orf | orf.track_ntseq_pos }.should == [9,18]
  end

  it "should combine a forward frame without ORFs in first seq" do
    s1 = "atggattaaatgta"
    #     ......---===xx
    s2 = "atggatttaattattataaa"
    #     x======xxx======xxx.
    fr = ShortFrameState.new s1,0,0
    fr.ntseq_pos.should == 0
    orfs = fr.get_stopstop_orfs
    orfs.size == 0 # in codons
    fr3 = FrameCodonHelpers::CreateShortFrame.create_right(fr,orfs,s2)
    fr3.ntseq_pos.should == 0
    fr3.codons.to_seq.should == "ATGGATTAAATGTAATGGATTTAATTATTATAA"
    norfs = fr3.get_stopstop_orfs
    orfs = norfs
    orfs.map{ | orf | orf.to_seq }.should == ["ATGTAA","TGGATTTAA","TTATTATAA"]
    orfs.map{ | orf | orf.track_ntseq_pos }.should == [9,9+6,9+6+9]
  end

  it "should combine a forward frame without ORFs in first seq" do
    s1 = "atggattaaatgta"
    #     ......---===xx
    s2 = "atggatttaatgtaaa"
    #     x======xxx
    fr = ShortFrameState.new s1,0,0
    fr.ntseq_pos.should == 0
    orfs = fr.get_stopstop_orfs
    orfs.size == 0 # in codons
    fr3 = FrameCodonHelpers::CreateShortFrame.create_right(fr,orfs,s2)
    # p fr3
    fr3.ntseq_pos.should == 0 # on the combined sequences
    fr3.codons.to_seq.should == "ATGGATTAAATGTAATGGATTTAATGTAAA"
    norfs = fr3.get_stopstop_orfs
    orfs += norfs
    orfs.map{ | orf | orf.to_seq }.should == ["ATGTAA", "TGGATTTAA"]
    orfs.map{ | orf | orf.track_ntseq_pos }.should == [9,9+6]
  end

  it "should combine a reverse frame" do
    s1 = "atggattaaatgta"
    s2 = "tatttaaatggatttaatgtaaatt"
    # now move the other way, as sequences get emitted on the left
    fr = ShortReversedFrameState.new s2,0,0
    # p fr
    fr.ntseq_pos.should == 0
    orfs = fr.get_stopstop_orfs
    # p orfs
    orfs.first.pos.should == 2 # in codons
    fr3 = FrameCodonHelpers::CreateShortFrame.create_right(fr,orfs,s2)
    fr3.ntseq_pos.should == 9
    fr3.codons.to_seq.should == "ATGTAATGGATTTAATGTAAA"
    norfs = fr3.get_stopstop_orfs
    orfs += norfs
    orfs.map{ | orf | orf.to_seq }.should == ["ATGGATTAA", "TGGATTTAA"]
    orfs.map{ | orf | orf.track_sequence_pos }.should == [0,11]
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
