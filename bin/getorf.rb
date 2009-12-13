#! /usr/bin/ruby
#
# Predict ORF's from nucleotide sequences using the EMBOSS translation routine
# and the BigBio predictors. The input is a fasta file, the output consists of
# a FASTA amino acid sequence file with matching nucleotide sequences
# (aa_heuristic.fa and nt_heuristic.fa respectively)
#
# You can choose the heuristic on the command line (default stopstop).
#
# (: pjotrp 2009 rblicense :)

USAGE =<<EOM
  ruby #{__FILE__} [-h startstop] inputfile(s)
EOM

$: << File.dirname(__FILE__)+'/../lib'

require 'biolib/emboss'
require 'bigbio'

if ARGV.size < 1
  print USAGE
  exit 1
end

heuristic = 'stopstop'

if ARGV[0] == '-h'
  heuristic = ARGV[1]
  ARGV.shift
  ARGV.shift
end

print "getorf #{heuristic}\n"
out = FastaPairedWriter.new('nt_'+heuristic+'.fa','aa_'+heuristic+'.fa')

ARGV.each do | fn |
  raise "File #{fn} does not exist" if !File.exist?(fn)
  nt = FastaReader.new(fn)
  trn_table = Biolib::Emboss.ajTrnNewI(1);

  nt.each do | rec |
    predict = PredictORF.new(rec.id,rec.descr,rec.seq,trn_table)
    orflist = predict.send(heuristic)
    orflist.each do | orf |
      out.write(orf.to_fastarec)
    end
  end
end





