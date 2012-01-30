
$: << "../lib"
ENV['DATA'] = '../test/data/EMBOSS'

require 'bigbio'

# Note that PredictORF, at this point, leaves trailing X's for the AA sequence

describe PredictORF, " when using a short simple nucleotide sequence" do
  before :all do 
    # initialize
    id = 'test'
    descr = 'Test'
    # sequence = 'AGCTGAATCTGGTAGATACCATCTTTAA'
    sequence = 'AGCTGAATCTGG'
    # trn_table = Biolib::Emboss.ajTrnNewI(1)
    trn_table = Bio::Big::TranslationAdapter.translation_table(1)
  
    @predictorf = PredictORF.new(id,descr,sequence,trn_table)
    @orflist = @predictorf.stopstop(0)
    # @orflist.each do | orf | p [orf.descr,orf] end
  end

  it "stopstop(0) should render six reading frames and seven ORF" do
    # >EMBOSS_001_1
    # S*IW
    # >EMBOSS_001_2
    # AESG
    # >EMBOSS_001_3
    # LNLX
    # >EMBOSS_001_4
    # PDSA
    # >EMBOSS_001_5
    # RFSX
    # >EMBOSS_001_6
    # QIQL
    @orflist[0].aa.seq.should == "S"
    @orflist[3].aa.seq[0..2].should == "LNL"
    @orflist[4].aa.seq[0..2].should == "PDS"
    @orflist.size.should == 7
  end

  # frame +1 - 4 codons S*IW 
  it "should give a first valid +1 frame" do
    orf = @orflist[1]
    orf.frame.should == 1
    orf.nt.start.should == 6
    orf.aa.seq.should == "IW"
    orf.nt.seq.should == "ATCTGG"
  end

  # frame +1 - 4 codons S*IW 
  it "should give a second valid +1 frame" do
    orf = @orflist[0]
    orf.frame.should equal 1
    orf.nt.start.should equal 0
    orf.aa.seq.should == "S"
    orf.nt.seq.should == "AGC"
  end

  # frame +2 - 3 codons AES
  it "should give a valid +2 frame" do 
    orf = @orflist[2]
    # pp @orflist
    # pp orf
    orf.nt.seq[0..8].should == "GCTGAATCT"
    orf.frame.should == 2
    orf.nt.start.should == 1
    # orf.nt.stop.should == 12 - EMBOS differs
    orf.aa.seq[0..2].should == "AES"
  end

  # frame +3 - 3 codons LNL
  it "should give a valid +3 frame" do
    orf = @orflist[3]
    orf.frame.should == 3
    orf.nt.start.should == 2
    # orf.nt.stop.should == 12
    orf.aa.seq[0..3].should == "LNL"
    orf.nt.seq[0..9].should == "CTGAATCTG"
  end

  # frame -1 - 4 codons PDSA 
  it "should give a valid -1 frame" do
    orf = @orflist[4]
    orf.frame.should == -1
    orf.nt.start.should == 0
    orf.nt.stop.should == 12
    orf.aa.seq.should == "PDSA"
    orf.nt.seq.should == "CCAGATTCAGCT"
  end

  # frame -2 - 3 codons RFSX
  it "should give a valid -3 frame" do
    orf = @orflist[5]
    orf.frame.should == -2
    orf.nt.start.should == 1
    # orf.nt.stop.should == 12
    orf.aa.seq[0..2].should == "RFS"
    orf.nt.seq[0..8].should == "CAGATTCAG"
  end

  # frame -3 - 3 codons QIQL
  it "should give a valid -3 frame" do
    orf = @orflist[6]
    orf.frame.should == -3
    orf.nt.start.should == 2
    # orf.nt.stop.should == 12
    orf.aa.seq[0..2].should == "QIQ"
    orf.nt.seq[0..8].should == "AGATTCAGC"
  end
end

describe PredictORF, " when using a more complicated nucleotide sequence" do
  before :all do 
    # initialize
    id = "PUT-157a-Arabidopsis_thaliana-126"
    descr = "PlantGDB Arabidopsis_thaliana Jan_15_2007"
    sequence = "ATCATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGT
                CTTCAGTGTACAGTATCAAGGGCTCGATCTGCGGTGGATGAGACATCAGATTCAGGAGCT
                TTTCAAAGAACTGCATCGACATCCGTAACTTCGTTTCAAAAGATTCCAATTCTCAGTTTC
                AGCTGAATCTGGTAGATACCATCTTTACATATCGTATGCTTGTCATGGGCTTCTAGATGC
                CTTTCATACTTAAAGATCAAAGGACTTGACGATGCAATAAGCTTCTCGTCTGTAAAACCC"
    # @trn_table = Biolib::Emboss.ajTrnNewI(1)
    @trn_table = Bio::Big::TranslationAdapter.translation_table(1)
    @predictorf = PredictORF.new(id,descr,sequence,@trn_table)
    @orflist = @predictorf.stopstop(0)
    # @orflist.each_with_index do | orf,i | p [i,orf.descr,orf.aa.seq,orf.nt.seq] end
    # >EMBOSS_001_1
    # IISNTSFLSLASKFTTRGSRLQCTVSRARSAVDETSDSGAFQRTASTSVTSFQKIPILSF
    # S*IW*IPSLHIVCLSWASRCLSYLKIKGLDDAISFSSVKP
    # >EMBOSS_001_2
    # SLATPASSLSLQSSLLVDLVFSVQYQGLDLRWMRHQIQELFKELHRHP*LRFKRFQFSVS
    # AESGRYHLYISYACHGLLDAFHT*RSKDLTMQ*ASRL*N
    # >EMBOSS_001_3
    # H*QHQLPLSRFKVHYSWISSSVYSIKGSICGG*DIRFRSFSKNCIDIRNFVSKDSNSQFQ
    # LNLVDTIFTYRMLVMGF*MPFILKDQRT*RCNKLLVCKT
    # >EMBOSS_001_4
    # GFYRREAYCIVKSFDL*V*KASRSP*QAYDM*RWYLPDSAETENWNLLKRSYGCRCSSLK
    # SS*I*CLIHRRSSP*YCTLKTRSTSSEL*SEREEAGVAND
    # >EMBOSS_001_5
    # VLQTRSLLHRQVL*SLSMKGI*KPMTSIRYVKMVSTRFS*N*ELESFETKLRMSMQFFEK
    # LLNLMSHPPQIEPLILYTEDEIHE**TLKRERGSWCC**X
    # >EMBOSS_001_6
    # GFTDEKLIASSSPLIFKYERHLEAHDKHTICKDGIYQIQLKLRIGIF*NEVTDVDAVL*K
    # APESDVSSTADRALDTVH*RRDPRVVNFEARERKLVLLMX
  end
  it "stopstop(0) should render ORFs" do
    @orflist[0].aa.seq[0..3].should == "IISN"
    @orflist[13].aa.seq[0..3].should == "GFYR"
    @orflist[22].aa.seq[0..3].should == "VLQT"
  end
  it "stopstop(0) should render 33 reading frames and seven ORF" do
    @orflist.size.should == 32
  end
  it "startstop(30) should render ORFs starting with a start codon" do
    orflist = @predictorf.startstop(5)
    # orflist.each do | orf | p [orf.descr,orf] end
    orflist[0].aa.seq.should == "MPFILKDQRT"
    orflist.size.should == 1
  end
  it "should never return an empty sequence" do
    orflist = @predictorf.stopstop(0)
    orflist.each do | orf |
      orf.nt.seq.size.should >= 0
    end
  end

  it "should return 3 sequences when the minsize is 132" do
    orflist = @predictorf.stopstop(44)
    orflist.size.should == 4
  end

  it "should return 2 sequences when the minsize is 133" do
    orflist = @predictorf.stopstop(45)
    orflist.size.should == 3
  end

  it "should have -1 frame" do
    sequence = "ATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGT"
    # >EMBOSS_001_4
    # TRSTSSEL*SEREEAGVAN
    predictorf = PredictORF.new('test','TEST',sequence,@trn_table)
    orflist = predictorf.stopstop(0)
    # orflist.each_with_index do | orf,i | p [i,orf.descr,orf.aa.seq,orf.nt.seq] end
    orflist[2].aa.seq[0..18].should == "QHQLPLSRFKVHYSWIS"
  end

  it "should correctly handle a sequence starting with a STOP codon" do
    sequence = "ATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGT"
    # >EMBOSS_001_3
    # *QHQLPLSRFKVHYSWIS
    predictorf = PredictORF.new('test','TEST',sequence,@trn_table)
    orflist = predictorf.stopstop(0)
    # orflist.each_with_index do | orf,i | p [i,orf.descr,orf.aa.seq,orf.nt.seq] end
    orflist[2].aa.seq[0..18].should == "QHQLPLSRFKVHYSWIS"
  end

end
