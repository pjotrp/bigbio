# ORF predictor class
#

require 'biolib/emboss'

# Helper class for storing ORF information
class ORFnucleotides
end

# Helper class for storing ORF information
class ORFaminoacids
end


class ORF
  attr_reader :nu, :aa
end

class PredictORF

  def initialize id, descr, seq, trn_table
    @id        = id
    @descr     = descr
    @seq       = seq
    @trn_table = trn_table
  end

  # Return a list of predicted ORFs with :minsize AA's
  def stopstop params = { :minsize => 30 }
    []
  end

end
