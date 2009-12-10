#! /usr/bin/ruby
#
# Translate nucleotide sequences into aminoacids sequences in all 
# reading frames.
#
#
# (: pjotrp 2009 rblicense :)

USAGE =<<EOM
  ruby #{__FILE__} inputfile(s)
EOM

$: << File.dirname(__FILE__)+'/../lib'

require 'biolib/emboss'
require 'bigbio'

if ARGV.size < 1
  print USAGE
  exit 1
end

ARGV.each do | fn |
  raise "File #{fn} does not exist" if !File.exist?(fn)
  nt = FastaReader.new(fn)
  trnTable = Biolib::Emboss.ajTrnNewI(1);

  nt.each { | rec |
      ajpseq   = Biolib::Emboss.ajSeqNewNameC(rec.seq,"Test sequence")

      [-3,-2,-1,1,2,3].each do | frame |
        ajpseqt  = Biolib::Emboss.ajTrnSeqOrig(trnTable,ajpseq,frame)
        aa       = Biolib::Emboss.ajSeqGetSeqCopyC(ajpseqt)
        print "> ",rec.id," ",frame.to_s,"\n"
        print aa,"\n"
    end
  }
end





