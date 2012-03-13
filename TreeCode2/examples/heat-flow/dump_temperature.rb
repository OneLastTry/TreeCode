#!/usr/bin/env ruby

require_relative '../3d-plasma-freq/parse_particles'
include Math

pos_opts = Options::parse
pos_opts.file = "positions.csv"
vel_opts = Options::parse
vel_opts.file = "velocities.csv"

positions = get_coords(pos_opts)
velocities = get_coords(vel_opts)

temperatures = {}
numbers = {}
bin_width = 0.25
positions.each_with_index do |pos,i|
	vel = velocities[i]
	bin = (pos.r / bin_width).floor
	temperatures[bin] = 0 if temperatures[bin] == nil
	numbers[bin] = 0 if numbers[bin] == nil
	temperatures[bin] += 0.5 * vel.r**2
	numbers[bin] += 1
end

#pp temperatures
#pp numbers

temperatures.each do |i,v|
	temperatures[i] /= numbers[i]
	volume = 4.0 /3 * PI * ( ((i+1)*bin_width)**3 - (i*bin_width)**3)
	density = numbers[i] / volume
	puts "#{i*bin_width}\t#{density}\t#{temperatures[i]}"
end
#pp temperatures
