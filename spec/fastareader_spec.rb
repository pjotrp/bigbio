
require 'rspec'

$: << "../lib"

require 'bigbio'

describe FastaGenomeReader, "when reading a full genome" do

  it "should load the genome file" do
    @genome = FastaGenomeReader.new("test/data/fasta/nt.fa", -> descr { return 'X',0,-1 }, 5)
    @genome[10].should == 'G'
    @genome.ref('X',112).should == 'A'
    @genome.ref('X',119).should == 'T'
    @genome.ref('X',120).should == 'C'
    @genome.ref('X',479).should == 'T'
    @genome.ref('X',480).should == 'T'
    @genome.ref('X',511).should == 'N'
    @genome.ref('X',560).should == 'T' # <- reads into the 3rd sequenc
  end

end
