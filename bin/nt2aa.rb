#! /usr/bin/ruby
#
# Translate nucleotide sequences into aminoacids sequences in all 
# reading frames.
#
USAGE =<<EOM
  ruby #{__FILE__} [--six-frame] inputfile(s)
EOM

$: << File.dirname(__FILE__)+'/../lib'

require 'bigbio'

if ARGV.size < 1
  print USAGE
  exit 1
end

do_sixframes = false
frames = [1]
if ARGV[0] == '--six-frame'
  ARGV.shift!
  do_sixframes = true
  frames = [-3,-2,-1,1,2,3]
end

require 'bigbio/adapters/translate'

ARGV.each do | fn |
  raise "File #{fn} does not exist" if !File.exist?(fn)
  nt = FastaReader.new(fn)
  trn_table = Bio::Big::TranslationAdapter.translation_table(1)
  
  nt.each { | rec |
      ajpseq   = Bio::Big::TranslationAdapter.pre_translate(rec.seq,"Test sequence")

      frames.each do | frame |
        aa  = Bio::Big::TranslationAdapter.translate(trn_table,frame, rec.seq, ajpseq)

        # ajpseqt  = Biolib::Emboss.ajTrnSeqOrig(trnTable,ajpseq,frame)
        # aa       = Biolib::Emboss.ajSeqGetSeqCopyC(ajpseqt)
        print "> ",rec.descr
        print " [",frame.to_s,"]" if do_sixframes
        print "\n"
        print aa,"\n"
    end
  }
end





