#!/usr/bin/env ruby

if ARGV.include? "-h"
	puts "Usage: #{__FILE__} dims num_circles radius num_electrons filename"
	exit
end

dims = ARGV[0].to_i
num_circles = ARGV[1].to_i
radius = ARGV[2].to_f
num_electrons = ARGV[3].to_i
file = ARGV[4]

puts "set size square"
puts "unset key"

str = "plot "
num_circles.times {|i| str << "'#{file}' u #{i*dims + 1}:#{i*dims + 2}:(#{radius}) every ::::1 w circles, "}
num_electrons.times {|i| str << "'#{file}' u #{i*dims + num_circles*2+1}:#{i*dims + num_circles*2+2}, "}
str.sub!(/, $/, '')

puts str


