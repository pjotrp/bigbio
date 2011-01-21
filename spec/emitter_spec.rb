
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
    fr = ShortFrameState.new s,minsize
    orfs = fr.get_stopstop_orfs
    orfs.map{ | orf | orf.to_seq }.should == []

    # Frame 1
    fr = ShortFrameState.new s[1..-1],minsize
    os = fr.get_stopstop_orfs
    os.map{ | orf | orf.to_seq }.should == ["TTNCGCGTGCCGCCTTCTTTCTCCTTTTTCTCTTTTACTTCTTCATCATCATCTTCTTCTTCTTCCTCTTCGATATTCGTCAGTGTGTGTATTTTGGGGAAAACTTTGTGA"]
    orfs += os
    # Frame 2
    fr = ShortFrameState.new s[2..-1],minsize
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
    fr = ShortFrameState.new s,minsize
    orfs = fr.get_startstop_orfs

    # Frame 1
    fr = ShortFrameState.new s[1..-1],minsize
    orfs += fr.get_startstop_orfs
    # Frame 2
    fr = ShortFrameState.new s[2..-1],minsize
    orfs += fr.get_startstop_orfs
    orfs.map{ | orf | orf.to_seq }.should == []
    orfs.size.should == 0
  end
end

describe Bio::Big::GetFrame, "when combining frames" do
  include Bio::Big
  it "should combine a forward frame" do
    seq_pos = 0
    s1 = "atggattaaatgta"
    s2 = "atggatttaatgtaaa"
    fr = ShortFrameState.new s1,0
    orfs = fr.get_stopstop_orfs
    orfs.last.rpos.should == 3
    remove = orfs.last.rpos*3
    remove.should == 9
    seq_pos += remove
    s3 = s1[remove..-1] + s2
    s3.should == "atgtaatggatttaatgtaaa"
    fr3 = fr.shortframe_right(s2)
    fr3 = ShortFrameState.new s3,0
    norfs = fr3.get_stopstop_orfs
    FrameCodonHelpers::TrackSequenceTrait.update_sequence_pos norfs, seq_pos

    orfs += norfs
    orfs.map{ | orf | orf.to_seq }.should == ["ATGGATTAA", "TGGATTTAA"]
    orfs.map{ | orf | orf.track_sequence_pos }.should == [nil,11]
    
  end

  it "should combine a reverse frame"
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
