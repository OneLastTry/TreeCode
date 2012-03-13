#!/usr/bin/env ruby

require 'pp'
require 'optparse'
require 'matrix'
require_relative 'parse_particles'


opts = Options::parse

get_coords(opts).each{|v| puts v.to_a.join("\t") }


