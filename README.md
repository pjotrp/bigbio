# BIGBIO

BigBio = BIG DATA for Ruby

BigBio is an initiative to a create high performance libraries for big data
computing in biology.

BigBio may use BioLib C/C++/D functions for increasing performance and
reducing memory consumption.

This is an experimental project. If you wish to contribute subscribe
to the BioRuby and/or BioLib mailing lists.

# Overview

* BigBio can translate nucleotide sequences to amino acid
  sequences using an EMBOSS C function, or BioRuby's translator.
* BigBio has an ORF emitter which parses DNA/RNA sequences and emits
  ORFs between START_STOP or STOP_STOP codons.
* BigBio has a FASTA file emitter, with iterates FASTA files and
  iterates sequences without loading everything in memory.

# Examples

## Iterate through a FASTA file

Read a file without loading the whole thing in memory

```ruby
fasta = FastaReader.new(fn)
fasta.each do | rec |
  print rec.descr,rec.seq
end
```

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

