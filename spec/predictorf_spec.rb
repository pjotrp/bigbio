
$: << "../lib"
ENV['DATA'] = '../test/data/EMBOSS'

require 'bigbio'

describe PredictORF, " when using a short simple nucleotide sequence" do
  before :all do 
    # initialize
    # id = "PUT-157a-Arabidopsis_thaliana-126"
    # descr = "PlantGDB Arabidopsis_thaliana Jan_15_2007"
    #  sequence = "ATCATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGT
    #             CTTCAGTGTACAGTATCAAGGGCTCGATCTGCGGTGGATGAGACATCAGATTCAGGAGCT
    #             TTTCAAAGAACTGCATCGACATCCGTAACTTCGTTTCAAAAGATTCCAATTCTCAGTTTC
    #             AGCTGAATCTGGTAGATACCATCTTTACATATCGTATGCTTGTCATGGGCTTCTAGATGC
    #             CTTTCATACTTAAAGATCAAAGGACTTGACGATGCAATAAGCTTCTCGTCTGTAAAACCC"
    id = 'test'
    descr = 'Test'
    # sequence = 'AGCTGAATCTGGTAGATACCATCTTTAA'
    sequence = 'AGCTGAATCTGG'
    trn_table = Biolib::Emboss.ajTrnNewI(1)
    @predictorf = PredictORF.new(id,descr,sequence,trn_table)
    @orflist = @predictorf.stopstop(0)
    @orflist.each do | orf | p [orf.descr,orf] end
  end

  it "stopstop(0) should render six reading frames and seven ORF" do
    @orflist = @predictorf.stopstop(0)
    @orflist.size.should == 7
  end

  # frame +1 - 4 codons S*IW 
  it "should give a first valid +1 frame" do
    orf = @orflist[5]
    orf.frame.should == 1
    orf.nt.start.should == 6
    orf.aa.seq.should == "IW"
    orf.nt.seq.should == "ATCTGG"
  end

  # frame +1 - 4 codons S*IW 
  it "should give a second valid +1 frame" do
    orf = @orflist[6]
    orf.frame.should equal 1
    orf.nt.start.should equal 0
    orf.aa.seq.should == "S"
    orf.nt.seq.should == "AGC"
  end

  # frame +2 - 3 codons AES
  it "should give a valid +2 frame" do 
    orf = @orflist[2]
    orf.frame.should == 2
    orf.nt.start.should == 1
    orf.nt.stop.should == 12
    orf.aa.seq.should == "AESG"
    orf.nt.seq.should == "GCTGAATCTGG"
  end

  # frame +3 - 3 codons LNL
  it "should give a valid +3 frame" do
    orf = @orflist[4]
    orf.frame.should == 3
    orf.nt.start.should == 2
    orf.nt.stop.should == 12
    orf.aa.seq.should == "LNLX"
    orf.nt.seq.should == "CTGAATCTGG"
  end

  # frame -1 - 4 codons PDSA 
  it "should give a valid -1 frame" do
    orf = @orflist[3]
    orf.frame.should == -1
    orf.nt.start.should == 0
    orf.nt.stop.should == 12
    orf.aa.seq.should == "PDSA"
    orf.nt.seq.should == "GGTCTAAGTCGA"
  end

  # frame -2 - 3 codons RFSX
  it "should give a valid -3 frame" do
    orf = @orflist[1]
    orf.frame.should == -2
    orf.nt.start.should == 1
    orf.nt.stop.should == 12
    orf.aa.seq.should == "RFSX"
    orf.nt.seq.should == "GTCTAAGTCGA"
  end

  # frame -3 - 3 codons QIQL
  it "should give a valid -2 frame" do
    orf = @orflist[0]
    orf.frame.should == -3
    orf.nt.start.should == 2
    orf.nt.stop.should == 12
    orf.aa.seq.should == "QIQL"
    orf.nt.seq.should == "TCTAAGTCGA"
  end
end
