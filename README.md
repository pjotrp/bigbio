# BIGBIO

BigBio = BIG DATA for Ruby

BigBio is an initiative to a create high performance libraries for big data
computing in biology.

BigBio may use BioLib C/C++/D functions for increasing performance and
reducing memory consumption.

In a way, this is an experimental project. I use it for
experimentation, but what is in here should work fine. If you wish to
contribute subscribe to the BioRuby and/or BioLib mailing lists
instead.

# Overview

* BigBio can translate nucleotide sequences to amino acid
  sequences using an EMBOSS C function, or BioRuby's translator.
* BigBio has an ORF emitter which parses DNA/RNA sequences and emits
  ORFs between START_STOP or STOP_STOP codons.
* BigBio has a FASTA file emitter, with iterates FASTA files and
  iterates sequences without loading everything in memory.

Warning: this software is experimental. Chech the issue list first.

# Examples

## Iterate through a FASTA file

Read a file without loading the whole thing in memory

```ruby
require 'bigbio'

fasta = FastaReader.new(fn)
fasta.each do | rec |
  print rec.descr,rec.seq
end
```

Since FastaReader parses the ID, write a tab file with id and sequence

```ruby
i = 1
print "num\tid\tseq\n"
FastaReader.new(fn).each do | rec |
  if rec.id =~ /(AT\w+)/
    print i,"\t",$1,"\t",rec.seq,"\n"
    i += 1
  end
end
```

wich, for example, can be turned into RDF with the
[bio-table](https://github.com/pjotrp/bioruby-table) biogem.

## Write a FASTA file

Write a FASTA file. The simple way

```ruby
fasta = FastaWriter.new(fn)
fasta.write('Test',"agtcta")
```

Any object can be passed in, however, as long as it responds to
'descr' and 'seq.to_s', or 'id' and 'seq.to_s'. E.g.

```ruby
class StorageObject
  attr_accessor :descr, :seq
end

mysequence = StorageObject.new
mysequence.descr = 'Test'
mysequence.seq = "agtcta"
```

write the FASTA file

```ruby
fasta = FastaWriter.new(fn)
fasta.write(mysequence)
```

## Transform a FASTA file

You can combine above FastaReader and FastaWriter to transform
sequences, e.g.

```ruby
fasta = FastaWriter.new(in_fn)
FastaReader.new(out_fn).each do | rec |
  # Strip the description down to the second ID
  (id1,id2) = /(\S+)\s+(\S+)/.match(rec.descr)
  fasta.write(id2,rec.seq)
end
```

The downside to this approach is the explicit file naming. What if you
want to use STDIN or some other source instead? I have come round to
the idea of using a combination of lambda and block. For example:

```ruby
  FastaReader::emit_fastarecord(-> {gets}) { |rec|
    print FastaWriter.to_fasta(rec)
  }
```

which takes STDIN line by line, and outputs FASTA on STDOUT. This is 
a correct design as the FastaReader and FastaWriter know nothing of
the mechanism fetching and displaying data. These can both be 'pure' 
functions. Note also that the data is never fully loaded into RAM.

Here the transformer functional style

```ruby
  FastaReader::emit_fastarecord(-> {gets}) { |rec|
    (id1,id2) = /(\S+)\s+(\S+)/.match(rec.descr)
    print FastaWriter.to_fasta(id2,req.seq)
  }
```

## Fetch ORFs from a sequence

BigBio can parse a sequence for ORFs. Together with the FastaReader
little memory gets used

```ruby
predictorf = PredictORF.new(id,descr,"ATCATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGT")
# get all ORFs between start and stop codons, longer than 30 bps
orfs = predictorf.startstop(30)
# get all sequences between stop codons
seqs = predictorf.stopstop(0)
```

## Rapid DNA/RNA to amino acid translation

Translate with EMBOSS C library, if linked, otherwise use BioRuby

```ruby
trn_table = Bio::Big::TranslationAdapter.translation_table(1)
translate = Nucleotide::Translate.new(trn_table)
aa_frames = translate.aa_6_frames("ATCATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGT")
```

# Command line tools

Some functionality comes also as executable command line tools (see the
./bin directory). Use the -h switch to get information. Current tools
are 

1. getorf: fetch all areas between start-stop and stop-stop codons in six frames (using EMBOSS when biolib is available)
2. nt2aa.rb: translate in six frames (using EMBOSS when biolib is available)

# Install

The easy way

```sh
gem install bio-bigbio
```

in your code

```ruby
require 'bigbio'
```

# Copyright

Copyright (c) 2011-2012 Pjotr Prins. See LICENSE for further details.

