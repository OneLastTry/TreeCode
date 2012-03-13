#!/usr/bin/env ruby

require 'pp'
require 'optparse'
require 'matrix'

class Options
	attr_accessor :dimensions, :timestep, :file

	def initialize
		@dimensions = 3
		@timestep = 0
		@file = nil
	end

	def self.parse
		return_opts = Options.new

		#Parse optional arguments
		opts = OptionParser.new do |opts|
			opts.banner = "Usage: #{$0} coord_file [OPTIONS]\n\t\"coord_file\" must be a csv file, every column a component of a particle's velocity or position."

			opts.on("-t", "--timetsep TIMESTEP", "Timestep at which to use results."){|t|
				return_opts.timestep = t.to_i
			}

			opts.on("-d", "--dimensions DIMS", "Number of dimensions we are working.."){|d|
				return_opts.dimensions = d.to_i
			}

	
			opts.on_tail("-h", "--help", "This help text."){|h|
				puts opts
				exit
			}
		end
		opts.parse!(ARGV)

		#Now, if we have more than one arg, failboat.
		#The remaining arg should be the database file
		if ARGV.length != 1
			STDERR.puts opts
			exit 1
		end
		return_opts.file = ARGV[0]
		return return_opts
	end
end

def get_coords(opts)
	coords = []
	File.open(opts.file, "r") do |f|
		opts.timestep.times {f.readline}
		cols = f.readline.split("\t")
		coord = []
		cols.each_with_index{|c,i|
			coord << c.to_f
			if (i+1) % opts.dimensions == 0
				v = Vector.elements(coord)
				if block_given?
					coords << v if yield(v)
				else
					coords << v
				end
				coord = []
			end
		}
	end
	return coords
end
