# Part of CHIMERA
# Messages container

class chimera_messages:
	def __init__(self):
		self.msg = {}
		self.log = []

		self.msg['sol/bounds'] = """Constructing solver with cylindric boundaries:
** left={0:.3g}, right={1:.3g}, radius={2:.3g}"""

		self.msg['sol/resols'] = """Spatial and temporal resolutions:
** dx={0:.3g}, dr={1:.3g}, dt={2:.3g}"""

		self.msg['sol/grdsz'] = """Grid sizes are:
** Nx={0:d}, Nr={1:d}, Mo={2:d}"""

		self.msg['sol/actv'] = """Solver is active from t={:g} to t={:g}"""

		self.msg['sol/rcut'] = """Rgrid is cut after {:.5g} (Nr is reduced to {:d})"""

		self.msg['sol/shift'] = """Spectral domain is shifted"""

		self.msg['sol/echo_intro'] = """Echo suppression is actived (expert feature)
** To change, set 'AntiEchoStrength' in solvers 'Features',
** or add 'NoAntiEcho' to solvers 'Features' to disactivate"""

		self.msg['sol/spcchrg'] = """Charge density will be considered"""

		self.msg['sol/stillspcs'] = """"Still" species will be treated as background"""

		self.msg['sol/pois'] = """Poisson correction will not be performed"""

	def print_(self, kw, args=None):
		if args is None:
			msg_str = self.msg[kw]
		else:
			msg_str =  self.msg[kw].format(*args)
		print(msg_str)
		self.log.append(msg_str)

msg = chimera_messages()
