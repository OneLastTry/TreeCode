#!/usr/bin/env ruby

include Math

def get_num
	while true
		u = rand()*2
		x = rand()
		return x if sin(2.0*PI*x)+1 > u
	end
end

100000.times {puts get_num }
