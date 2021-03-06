# splotter
Python based S-Parameter Plotting Tool
authors: Michael Cracraft, Andrew Becker (QT Updates)

Usage: (No command line arguments exist at this time.)
splotter (or Splotter_PyQt_lnx.py)

In *NIX OSes, splotter is included as a symbolic link to Splotter_PyQt_lnx.py.

Notes: 
The JSON replotting feature from files appears to be broken in case anyone feels
like fixing it before I get time to do so.

Prerequisites: 
splotter is currently working well with an Anaconda installation of Python 3.5
on Linux.  It is not thoroughly tested on Windows, but it should function if the
following prereqs are present in addition to the modules included in the repo.
- numpy
- matplotlib
- qt (Qt5)

Below is the response to "conda list" for my Python 3.5 environment as of 6
January 2017.  This constitutes a sufficient, but not necessary list of
packages.

# packages in environment at /home/xxxxxxxx/anaconda3:
#
_license                  1.1                      py35_1  
_nb_ext_conf              0.3.0                    py35_0  
alabaster                 0.7.9                    py35_0  
anaconda                  custom                   py35_0  
anaconda-clean            1.1.0                    py35_0  
anaconda-client           1.6.0                    py35_0  
anaconda-navigator        1.3.2                    py35_0  
argcomplete               1.0.0                    py35_1  
astroid                   1.4.7                    py35_0  
astropy                   1.3                 np111py35_0  
babel                     2.3.4                    py35_0  
backports                 1.0                      py35_0  
beautifulsoup4            4.5.3                    py35_0  
bitarray                  0.8.1                    py35_0  
blaze                     0.10.1                   py35_0  
bokeh                     0.12.3                   py35_0  
boto                      2.45.0                   py35_0  
bottleneck                1.1.0               np111py35_0  
cairo                     1.14.6                        0  
cffi                      1.9.1                    py35_0  
chest                     0.2.3                    py35_0  
click                     6.6                      py35_0  
cloudpickle               0.2.1                    py35_0  
clyent                    1.2.2                    py35_0  
colorama                  0.3.7                    py35_0  
conda                     4.2.13                   py35_0  
conda-build               2.1.0                    py35_0  
conda-env                 2.6.0                         0  
conda-verify              2.0.0                    py35_0  
configobj                 5.0.6                    py35_0  
contextlib2               0.5.4                    py35_0  
cryptography              1.7.1                    py35_0  
curl                      7.49.0                        1  
cycler                    0.10.0                   py35_0  
cython                    0.25.2                   py35_0  
cytoolz                   0.8.2                    py35_0  
dask                      0.12.0                   py35_0  
datashape                 0.5.4                    py35_0  
dbus                      1.10.10                       0  
decorator                 4.0.10                   py35_1  
dill                      0.2.5                    py35_0  
docutils                  0.12                     py35_2  
dynd-python               0.7.2                    py35_0  
entrypoints               0.2.2                    py35_0  
et_xmlfile                1.0.1                    py35_0  
expat                     2.1.0                         0  
fastcache                 1.0.2                    py35_1  
filelock                  2.0.7                    py35_0  
flask                     0.12                     py35_0  
flask-cors                2.1.2                    py35_0  
fontconfig                2.12.1                        0  
freetype                  2.5.5                         1  
geopy                     1.11.0                    <pip>
get_terminal_size         1.0.0                    py35_0  
gevent                    1.2.0                    py35_0  
glib                      2.50.2                        0  
greenlet                  0.4.11                   py35_0  
gst-plugins-base          1.8.0                         0  
gstreamer                 1.8.0                         0  
h5py                      2.6.0               np111py35_2  
harfbuzz                  0.9.39                        1  
hdf5                      1.8.17                        1  
heapdict                  1.0.0                    py35_1  
ibm-db                    2.0.7                     <pip>
icu                       54.1                          0  
idna                      2.1                      py35_0  
imagesize                 0.7.1                    py35_0  
ipykernel                 4.5.2                    py35_0  
ipython                   5.1.0                    py35_0  
ipython_genutils          0.1.0                    py35_0  
ipywidgets                5.2.2                    py35_0  
itsdangerous              0.24                     py35_0  
jbig                      2.1                           0  
jdcal                     1.3                      py35_0  
jedi                      0.9.0                    py35_1  
jinja2                    2.8                      py35_1  
jpeg                      8d                            2  
jsonschema                2.5.1                    py35_0  
jupyter                   1.0.0                    py35_3  
jupyter_client            4.4.0                    py35_0  
jupyter_console           5.0.0                    py35_0  
jupyter_core              4.2.1                    py35_0  
lazy-object-proxy         1.2.1                    py35_0  
libdynd                   0.7.2                         0  
libffi                    3.2.1                         1  
libgcc                    5.2.0                         0  
libgfortran               3.0.0                         1  
libiconv                  1.14                          0  
libpng                    1.6.22                        0  
libsodium                 1.0.10                        0  
libtiff                   4.0.6                         2  
libxcb                    1.12                          1  
libxml2                   2.9.4                         0  
libxslt                   1.1.29                        0  
llvmlite                  0.15.0                   py35_0  
locket                    0.2.0                    py35_1  
lxml                      3.7.0                    py35_0  
markupsafe                0.23                     py35_2  
matplotlib                1.5.3               np111py35_1  
mistune                   0.7.3                    py35_0  
mkl                       11.3.3                        0  
mkl-service               1.1.2                    py35_2  
mpmath                    0.19                     py35_1  
multipledispatch          0.4.9                    py35_0  
nb_anacondacloud          1.2.0                    py35_0  
nb_conda                  2.0.0                    py35_0  
nb_conda_kernels          2.0.0                    py35_0  
nbconvert                 4.2.0                    py35_0  
nbformat                  4.2.0                    py35_0  
nbpresent                 3.0.2                    py35_0  
networkx                  1.11                     py35_0  
nltk                      3.2.1                    py35_0  
nose                      1.3.7                    py35_1  
notebook                  4.3.0                    py35_0  
numba                     0.30.0              np111py35_0  
numexpr                   2.6.1               np111py35_1  
numpy                     1.11.2                   py35_0  
odo                       0.5.0                    py35_1  
openpyxl                  2.4.0                    py35_0  
openssl                   1.0.2j                        0  
pandas                    0.19.2              np111py35_0  
partd                     0.3.6                    py35_0  
patchelf                  0.9                           0  
path.py                   9.0.1                    py35_0  
pathlib2                  2.1.0                    py35_0  
patsy                     0.4.1                    py35_0  
pcre                      8.39                          1  
pep8                      1.7.0                    py35_0  
pexpect                   4.0.1                    py35_0  
pickleshare               0.7.4                    py35_0  
pillow                    3.4.2                    py35_0  
pip                       9.0.1                    py35_1  
pip                       9.0.1                     <pip>
pixman                    0.34.0                        0  
pkginfo                   1.4.1                    py35_0  
ply                       3.9                      py35_0  
prompt_toolkit            1.0.9                    py35_0  
psutil                    5.0.1                    py35_0  
ptyprocess                0.5.1                    py35_0  
py                        1.4.31                   py35_0  
pyasn1                    0.1.9                    py35_0  
pycosat                   0.6.1                    py35_1  
pycparser                 2.17                     py35_0  
pycrypto                  2.6.1                    py35_4  
pycurl                    7.43.0                   py35_0  
pyflakes                  1.4.0                    py35_0  
pygments                  2.1.3                    py35_0  
pylint                    1.5.4                    py35_1  
pyodbc                    3.0.10                   py35_1  
pyopenssl                 16.2.0                   py35_0  
pyparsing                 2.1.4                    py35_0  
pyqt                      5.6.0                    py35_1  
pytables                  3.3.0               np111py35_0  
pytest                    3.0.5                    py35_0  
python                    3.5.2                         0  
python-dateutil           2.6.0                    py35_0  
pytz                      2016.10                  py35_0  
PyVISA                    1.8                       <pip>
pyyaml                    3.12                     py35_0  
pyzmq                     16.0.2                   py35_0  
qt                        5.6.2                         2  
qtawesome                 0.4.1                    py35_0  
qtconsole                 4.2.1                    py35_1  
qtpy                      1.1.2                    py35_0  
readline                  6.2                           2  
redis                     3.2.0                         0  
redis-py                  2.10.5                   py35_0  
requests                  2.12.4                   py35_0  
rope                      0.9.4                    py35_1  
ruamel_yaml               0.11.14                  py35_0  
scikit-image              0.12.3              np111py35_1  
scikit-learn              0.18.1              np111py35_0  
scipy                     0.18.1              np111py35_0  
setuptools                27.2.0                   py35_0  
simplegeneric             0.8.1                    py35_1  
singledispatch            3.4.0.3                  py35_0  
sip                       4.18                     py35_0  
six                       1.10.0                   py35_0  
snowballstemmer           1.2.1                    py35_0  
sockjs-tornado            1.0.3                    py35_0  
sphinx                    1.5.1                    py35_0  
sphinx_rtd_theme          0.1.9                    py35_0  
spyder                    3.0.2                    py35_0  
sqlalchemy                1.1.4                    py35_0  
sqlite                    3.13.0                        0  
statsmodels               0.6.1               np111py35_1  
sympy                     1.0                      py35_0  
terminado                 0.6                      py35_0  
tk                        8.5.18                        0  
toolz                     0.8.2                    py35_0  
tornado                   4.4.2                    py35_0  
traitlets                 4.3.1                    py35_0  
unicodecsv                0.14.1                   py35_0  
unixodbc                  2.3.4                         0  
wcwidth                   0.1.7                    py35_0  
werkzeug                  0.11.13                  py35_0  
wheel                     0.29.0                   py35_0  
widgetsnbextension        1.2.6                    py35_0  
wrapt                     1.10.8                   py35_0  
xlrd                      1.0.0                    py35_0  
xlsxwriter                0.9.6                    py35_0  
xlwt                      1.1.2                    py35_0  
xz                        5.2.2                         1  
yaml                      0.1.6                         0  
zeromq                    4.1.5                         0  
zlib                      1.2.8                         3  
