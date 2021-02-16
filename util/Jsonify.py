import json

class JsonifySession(object):

	def __init__(self, session_name='Session.json'):
		self.session = session_name
		self.plot_dictionary = dict()
		self.plot_dictionary['num_plots'] = 0

	def add_plot(self, plot_dict):
		"""
		Adds a plot to the session dictionary
		"""
		num = self.plot_dictionary['num_plots']
		label = 'plot_%d' % (num,)
		self.plot_dictionary['num_plots'] = num + 1
		self.plot_dictionary[label] = plot_dict

	def get_session(self):
		"""
		Returns the session dictionary as jsonified text
		"""
		return json.dumps(self.plot_dictionary)

	def save(self):
		"""
		Writes the Jsonified session to a file, Specified by the session name.
		"""
		json_text = self.get_session()
		with open(self.session, 'w') as json_file:
			json_file.write(json_text)


class PlotDict(dict):

	def __init__(self, title, ylabel, xlabel):
		super(PlotDict, self).__init__()
		self['title'] = title
		self['ylabel'] = ylabel
		self['xlabel'] = xlabel
		self['num_lines'] = 0

	def add_line(self, xarray, yarray, label):
		num = self['num_lines']
		line_label = 'line_%d' % (num,)
		self['num_lines'] = num + 1
		line_dict = {
			'x': xarray, 
			'y': yarray,
			'label': label
		}
		self[line_label] = line_dict

	def set_limits(self, x_lim, y_lim):
		self['x_lim'] = x_lim
		self['y_lim'] = y_lim
		

if __name__ == '__main__':
	"""
	This is just simple testing of the module
	"""
	that = PlotDict('plot that', 'y', 'x')
	this = PlotDict('plot this', 'y', 'x')
	that.add_line([1,2,3], [1,2,3], 'line1')
	that.add_line([1,2,3], [3,2,1], 'line2')
	this.add_line([1,2,3], [3,2,1], 'only line')
	session = JsonifySession()
	session.add_plot(that)
	session.add_plot(this)
	session.save()
