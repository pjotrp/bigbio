
require 'singleton'

module Bio
  module Big
    class Environment
      include Singleton

      attr_accessor :log, :biolib
    end
  end
end
