class Interval:
	
	include_strings = ['i', 'inclusive', 'include']
	exclude_strings = ['e', 'exclusive', 'exclude']
	valid_type_strings = include_strings + exclude_strings
	
	def __init__(self, *args, **kwargs):
		
		if len(args) == 1:
			lower, upper = args[0].split(',') # Should error out if bad formatting
			self.lower_bound_type = 'e' if lower[0] == '(' else 'i'
			self.lower_bound = float(lower.strip()[1:])
			self.upper_bound_type = 'e' if upper[-1] == ')' else 'i'
			self.upper_bound = float(upper.strip()[:-1])
		else:
			# Otherwise, assume 4 arguments. This should error out if there
			# are an incorrect number of arguments. Tacky but yea
			lbound_type, lbound, ubound, ubound_type = args
			
			self.lower_bound = lbound
			self.upper_bound = ubound
			
			if lbound_type in self.valid_type_strings:
				self.lower_bound_type = lbound_type
			else:
				print("Invalid bound type! Defaulting to inclusive")
				self.lower_bound_type = 'inclusive'
			
			if ubound_type in self.valid_type_strings:
				self.upper_bound_type = ubound_type
			else:
				print("Invalid bound type! Defaulting to inclusive")
				self.upper_bound_type = 'inclusive'
	
	def contains(self, val):
		if (self.lower_bound_type in self.include_strings):
			lower_good = (val >= self.lower_bound)
		else:
			lower_good = (val > self.lower_bound)
		
		if (self.upper_bound_type in self.include_strings):
			upper_good = (val <= self.upper_bound)
		else:
			upper_good = (val < self.upper_bound)
		
		return (lower_good and upper_good)