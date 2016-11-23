

import os
import sys
ROOTPATH = os.path.dirname(os.path.dirname(__file__,))

if sys.platform == 'win32': #Windows
    BABEL_COM = 'C:\\"Program Files (x86)"\\OpenBabel-2.3.2\\babel '
    DELETE_COM = 'del '
    PYTHON_COM = 'python '
elif sys.platform == 'darwin': #OS X
	BABEL_COM = 'babel'
	DELETE_COM = 'rm '
	PYTHON_COM = 'python '
else: #Linux
    #BABEL_COM = '/usr/local/bin/babel '
    BABEL_COM = '~/Tools/openbabel-install/bin/babel  '
    DELETE_COM = 'rm '
    PYTHON_COM = 'python '