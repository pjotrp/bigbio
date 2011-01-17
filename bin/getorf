#! /usr/bin/ruby
#
# Predict ORF's from nucleotide sequences using the BigBio predictors. 
# The input is a fasta file, the output consists of
# a FASTA amino acid sequence file with matching nucleotide sequences
# (aa_heuristic.fa and nt_heuristic.fa respectively)
#
# You can choose the heuristic on the command line (default stopstop).
#
# Author:: Pjotr Prins
# Copyright:: 2009-2011
# License:: Ruby License
#
# Copyright (C) 2009-2011 Pjotr Prins <pjotr.prins@thebird.nl>

rootpath = File.dirname(File.dirname(__FILE__))
$: << File.join(rootpath,'lib')

BIGBIO_VERSION = File.new(File.join(rootpath,'VERSION')).read.chomp

USAGE =<<EOM
  ruby #{__FILE__} [-h stopstop] [--min-size 30] inputfile(s)
EOM

# require 'biolib/emboss'
require 'bigbio'
require 'optparse'

$stderr.print "getorf BioRuby BigBio Plugin "+BIGBIO_VERSION+" Copyright (C) 2009-2011 Pjotr Prins <pjotr.prins@thebird.nl>\n\n"

heuristic = 'stopstop'
minsize   = 30
longest_match = false

opts = OptionParser.new() { |opts|
  opts.on_tail("-?", "--help", "Print this message") {
    print(USAGE)
    print(opts)
    print <<EXAMPLE
   
EXAMPLE
    exit()
  }
   
  opts.on("-h heuristic", String, "Heuristic (stopstop)") do | s |
    heuristic = s
  end
  opts.on("-s size", "--min-size", Integer, "Minimal sequence size") do | n |
    minsize = n
  end
  opts.on("--longest", "Only get longest ORF match") do 
    longest_match = true
  end
}
opts.parse!(ARGV)
if ARGV.size == 0
  print USAGE
  exit 1
end

print "Heuristic is #{heuristic}\n"
print "Minsize #{minsize}\n"

out = FastaPairedWriter.new('nt_'+heuristic+'.fa','aa_'+heuristic+'.fa')

ARGV.each do | fn |
  raise "File #{fn} does not exist" if !File.exist?(fn)
  nt = FastaReader.new(fn)
  trn_table = Biolib::Emboss.ajTrnNewI(1);

  nt.each do | rec |
    predict = PredictORF.new(rec.id,rec.descr,rec.seq,trn_table)
    if longest_match
      orf = predict.send('longest_'+heuristic,minsize)
      out.write(orf.to_fastarec) if orf
    else
      orflist = predict.send(heuristic,minsize)
      orflist.each do | orf |
        out.write(orf.to_fastarec)
      end
    end
  end
end





