import json
import sys
import matplotlib as mpl
mpl.use('qt4agg')
import matplotlib.pyplot as plt


def load_json_file(path):
	with open(path, 'r') as file:
		json_text = file.read()
	return json_text


def plot_json(json_session):
	session_dict = json.loads(json_session)
	num_plots = session_dict['num_plots']
	for i in range(num_plots):
		fig = plt.figure()
		axes = fig.add_subplot(1,1,1)
		label = 'plot_%d' % (i,)
		plot_dict = session_dict[label]
		num_lines = plot_dict['num_lines']
		for j in range(num_lines):
			label_ = 'line_%d' % (j,)
			line = plot_dict[label_]
			x = line['x']
			y = line['y']
			axes.plot(x, y, label=line['label'], picker=5)
		axes.grid(True)
		axes.legend(loc='best')
		axes.set_xlabel(plot_dict['xlabel'])
		axes.set_ylabel(plot_dict['ylabel'])
		axes.set_title(plot_dict['title'])
		axes.set_xlim(plot_dict['x_lim'])
		axes.set_ylim(plot_dict['y_lim'])
		fig.show()  # Keeps Each plot in it's own window
	plt.show()  # Blocks so the plots stay up until the user closes them


if __name__ == '__main__':
	argv = sys.argv[1:]
	if len(argv) < 1:
		print( 'Please enter path to the json file you would like to plot.' )
	elif len(argv) > 1:
		print( 'Too many command line arguments. Please only enter the path to the json flie.' )
	else:
		path = argv[0]
		text = load_json_file(path)
		plot_json(text)



