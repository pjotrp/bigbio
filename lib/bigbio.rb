# BigBio libraries

require 'bigbio/environment'

# find local plugin installation, and use it when there
rootpath = File.dirname(File.dirname(__FILE__))
bio_logger_path = File.join(rootpath,'..','bioruby-logger','lib')
if File.directory? bio_logger_path
  $: << bio_logger_path
  $stderr.print "bio-logger loaded directly\n"
else
  require "rubygems"
  gem "bio-logger"
end
require 'bio-logger'

log = Bio::Log::LoggerPlus.new('bigbio')
Bio::Big::Environment.instance.log = log

begin
  require 'biolib/emboss'
  Bio::Big::Environment.instance.biolib = true
rescue LoadError
  log.outputters = Bio::Log::Outputter.stderr
  log.warn "BioLib functionality not loaded"
end


# module BigBio

  autoload :FastaReader, 'bigbio/db/fasta'
  autoload :FastaWriter, 'bigbio/db/fasta'
  autoload :FastaPairedReader, 'bigbio/db/fasta'
  autoload :FastaPairedWriter, 'bigbio/db/fasta'
  autoload :BlastClust, 'bigbio/db/blast'
  autoload :PredictORF, 'bigbio/sequence/predictorf'

# end
