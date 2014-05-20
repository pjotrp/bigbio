
require 'rspec'

$: << "../lib"

require 'bigbio'

describe FastaGenomeReader, "when reading a full genome" do

  it "should load the genome file" do
    FastaGenomeReader.new("test/data/fasta/nt.fa", -> {})
  end

end
