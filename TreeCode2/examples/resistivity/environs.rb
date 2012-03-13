#!/usr/bin/env ruby
include Math
require 'irb'

M_e = 9.10938188e-31
E_0 = 8.85418782e-12
E   = 1.60217646e-19

def t_c(r_c)
	return sqrt(4.0 * PI * M_e * E_0 * r_c**3 / E**2)
end

ARGV.clear # otherwise all script parameters get passed to IRB
IRB.start
