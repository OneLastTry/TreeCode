#!/usr/bin/env ruby

require 'pp'
require 'optparse'
require 'matrix'
require_relative 'parse_particles'


opts = Options::parse

(0 .. 10000).step(1) do |t|
	opts.timestep = t
	coords = get_coords(opts) {|c| c[0] > -0.1 && c[0] < 0.1}
	#coords.each{|v| puts v.to_a.join("\t") }
	puts "#{t}\t#{coords.length}"
end

